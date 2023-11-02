#pragma once
#include <random>
#include "upcxx/upcxx.hpp"

constexpr char distribution = 0; // 0 -> Cyclic, 1 -> Contiguous

class Tile
{
public:
    Tile(uint64_t row_min, uint64_t col_min, uint64_t row_max, uint64_t col_max)
    : row_min(row_min), col_min(col_min), row_max(row_max), col_max(col_max) {}

    bool operator<(const Tile& other) const noexcept
    {
        return row_max < other.row_min ? true : row_min > other.row_max ? false : col_max < other.col_min ? true : col_min > other.col_max ? false : false;
    }

    uint64_t row_min, col_min, row_max, col_max;
};

class Darray // 2D
{
public:
    using value_type = double;

    Darray( const upcxx::team& t = upcxx::world() ) : team(t) {} // TODO: default(ish), if we want tp preclude this behavior, then we have to add some extra logic to the drivers  
    Darray(const std::vector<uint64_t>& tile_heights, const std::vector<uint64_t>& tile_widths, const upcxx::team& team = upcxx::world()) noexcept :
        tile_heights(tile_heights), tile_widths(tile_widths), team(team)
        {
            nrows = std::accumulate(tile_heights.begin(), tile_heights.end(), 0);
            ncols = std::accumulate(tile_widths.begin(), tile_widths.end(), 0);
            nelems = nrows * ncols;
            ntiles = tile_heights.size() * tile_widths.size();
            local_nelems = 0;
        }

    void allocate() noexcept
    {
        auto ranks = team.rank_n();
        uint64_t tiles_per_proc = (ntiles + ranks - 1) / ranks;
        upcxx::intrank_t owner = 0;
        std::vector<uint64_t> offsets(ranks, 0);

        for(uint64_t i = 0, ii = 0; i < nrows; i += tile_heights[ii], ++ii)
            for(uint64_t j = 0, jj = 0; j < ncols; j += tile_widths[jj], ++jj)
            {
                tiles[{i, j, i + tile_heights[ii] - 1, j + tile_widths[jj] - 1}] = {owner, offsets[owner]};

                offsets[owner] += tile_heights[ii] * tile_widths[jj];

                if constexpr(distribution == 0) // Cyclic
                    owner = (owner + 1) % ranks;
                else
                    owner = tiles.size() / tiles_per_proc;
            }

        local_nelems = offsets[team.rank_me()];

        local_gptr = upcxx::new_array<value_type>(local_nelems);

        gptrs.resize(ranks);
        upcxx::promise<> p(ranks);
        for(auto r = 0; r < ranks; ++r)
            upcxx::broadcast(local_gptr, r, team).then([this, &p, r](auto result)
            {
                gptrs[r] = result;
                p.fulfill_anonymous(1);
            });
        p.get_future().wait();
    }

    // Get a chunk of contiguous elements. e.g. a single element, a whole tile, or a contiguous
    // subset of a tile.
    template <bool copy = true>
    void get_contig(uint64_t row_min, uint64_t col_min, uint64_t row_max, uint64_t col_max, value_type** buf) const noexcept
    {
        const auto& t = tiles.find({row_min, col_min, row_max, col_max});

        if(t != tiles.end())
        {
            auto offset = (row_min - t->first.row_min) * (t->first.col_max - t->first.col_min + 1) + (col_min - t->first.col_min);
            auto remote_addr = gptrs[t->second.first] + t->second.second + offset;
            auto count = (row_max - row_min + 1) * (col_max - col_min + 1);

            if(remote_addr.is_local())
                if constexpr(copy)
                    std::copy_n(remote_addr.local(), count, *buf);
                else
                    *buf = remote_addr.local();
            else
                upcxx::rget(remote_addr, *buf, count).wait();
        }
    }

    // Put a chunk of contiguous elements. e.g. a single element, a whole tile, or a contiguous
    // subset of a tile.
    void put_contig(uint64_t row_min, uint64_t col_min, uint64_t row_max, uint64_t col_max, const value_type* buf) const noexcept
    {
        const auto& t = tiles.find({row_min, col_min, row_max, col_max});

        if(t != tiles.end())
        {
            auto offset = (row_min - t->first.row_min) * (t->first.col_max - t->first.col_min + 1) + (col_min - t->first.col_min);
            auto remote_addr = gptrs[t->second.first] + t->second.second + offset;
            auto count = (row_max - row_min + 1) * (col_max - col_min + 1);

            if(remote_addr.is_local())
                std::copy_n(buf, count, remote_addr.local());
            else
                upcxx::rput(buf, remote_addr, count).wait();
        }
    }

    // TODO: This function will be improved by using logic + rget_strided(...)
    void get(uint64_t row_min, uint64_t col_min, uint64_t row_max, uint64_t col_max, value_type* buf) const noexcept
    {
        uint64_t next = 0;
        upcxx::promise<> p;
        std::unordered_map<upcxx::intrank_t, std::pair<std::vector<upcxx::global_ptr<value_type>>, std::vector<value_type*>>> all_gets;

        for(uint64_t i = row_min; i <= row_max; ++i)
            for(uint64_t j = col_min; j <= col_max; ++j)
            {
                const auto& t = tiles.find({i, j, i, j});
                auto offset = (i - t->first.row_min) * (t->first.col_max - t->first.col_min + 1) + (j - t->first.col_min);
                auto remote_addr = gptrs[t->second.first] + t->second.second + offset;
                auto local_addr = buf + next++;

                if(remote_addr.is_local())
                    *local_addr = *remote_addr.local();
                else {
                    all_gets[t->second.first].first.push_back(remote_addr);
                    all_gets[t->second.first].second.push_back(local_addr);
                }
            }

        auto sz = sizeof(value_type);

        for(const auto& x: all_gets)
            upcxx::rget_regular(x.second.first.begin(), x.second.first.end(), sz,
                                x.second.second.begin(), x.second.second.end(), sz,
                                upcxx::operation_cx::as_promise(p));

        p.finalize().wait();
    }

    // TODO: This function will be improved by using logic + rput_strided(...)
    void put(uint64_t row_min, uint64_t col_min, uint64_t row_max, uint64_t col_max, value_type* buf) const noexcept
    {
        uint64_t next = 0;
        std::unordered_map<upcxx::intrank_t, std::pair<std::vector<value_type*>, std::vector<upcxx::global_ptr<value_type>>>> all_puts;

        for(uint64_t i = row_min; i <= row_max; ++i)
            for(uint64_t j = col_min; j <= col_max; ++j)
            {
                const auto& t = tiles.find({i, j, i, j});
                auto offset = (i - t->first.row_min) * (t->first.col_max - t->first.col_min + 1) + (j - t->first.col_min);
                auto remote_addr = gptrs[t->second.first] + t->second.second + offset;
                auto local_addr = buf + next++;

                if(remote_addr.is_local())
                    *remote_addr.local() = *local_addr;
                else {
                    all_puts[t->second.first].first.push_back(local_addr);
                    all_puts[t->second.first].second.push_back(remote_addr);
                }
            }

        upcxx::promise<> p;
        auto sz = sizeof(value_type);

        for(const auto& x: all_puts)
        upcxx::rput_regular(x.second.first.begin(), x.second.first.end(), sz,
                            x.second.second.begin(), x.second.second.end(), sz,
                            upcxx::operation_cx::as_promise(p));

        p.finalize().wait();
    }

    void fill() noexcept
    {
        std::mt19937 gen(team.rank_me());
        std::uniform_real_distribution<value_type> distribution(0.0, 1.0);
        std::generate_n(local_gptr.local(), local_nelems, [&](){ return distribution(gen); });
        //std::fill_n(local_gptr.local(), local_nelems, team.rank_me());
    }

    void print_local() const noexcept
    {
        auto ptr = local_gptr.local();
        for(uint64_t i = 0; i < local_nelems; ++i)
        {
            std::cout << std::fixed << *ptr << ", ";
            ++ptr;
        }
        std::cout << std::endl;
    }

    void print(std::ostream& out = std::cout) const noexcept
    {
        if(!team.rank_me())
        {
            std::vector<value_type> values;
            // value_type* values;

            for(const auto& t : tiles)
            {
                values.resize((t.first.row_max - t.first.row_min + 1) * (t.first.col_max - t.first.col_min + 1));

                get_contig(t.first.row_min, t.first.col_min, t.first.row_max, t.first.col_max, (value_type**)&values);
                // OR:
                // get_contig<false>(t.first.row_min, t.first.col_min, t.first.row_max, t.first.col_max, &values);

                uint64_t k = 0;

                for(uint64_t i = t.first.row_min; i <= t.first.row_max; ++i)
                    for(uint64_t j = t.first.col_min; j <= t.first.col_max; ++j)
                        out << '(' << i << ',' << j << ") -> " << std::fixed << values[k++] << std::endl;
            }
        }
    }

    template <char operation = '='>
    Darray& op(const value_type* buf) noexcept
    {
        for(const auto& t : tiles)
            if(t.second.first == team.rank_me())
            {
                auto ptr = local_gptr.local() + t.second.second;
                auto t_width = t.first.col_max - t.first.col_min + 1;
                auto t_size = (t.first.row_max - t.first.row_min + 1) * t_width;
                auto global_idx = t.first.row_min * ncols + t.first.col_min;
                auto jump = ncols - t_width;

                for(uint64_t i = 1; i <= t_size; ++i)
                {
                    if constexpr(operation == '*')
                        *ptr++ *= *(buf + global_idx++);
                    else if constexpr(operation == '+')
                        *ptr++ += *(buf + global_idx++);
                    else if constexpr(operation == '=')
                        *ptr++ = *(buf + global_idx++);
                    if(i % t_width == 0)
                        global_idx += jump;
                }
            }

        return *this;
    }

    void deallocate() noexcept
    {
        local_nelems = 0;
        tiles.clear();
        gptrs.clear();
        upcxx::barrier(team);
        upcxx::delete_array(local_gptr);
        upcxx::barrier(team);
    }

    uint64_t nrows, ncols, nelems, ntiles, local_nelems;
    std::vector<uint64_t> tile_heights, tile_widths;
    std::vector<upcxx::global_ptr<value_type>> gptrs;
    std::map<Tile, std::pair<upcxx::intrank_t, uint64_t>> tiles;
    upcxx::global_ptr<value_type> local_gptr;
    const upcxx::team& team;
};
