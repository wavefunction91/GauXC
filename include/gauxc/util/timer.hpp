#pragma once
#include <chrono>
#include <map>
#include <string>
#include <type_traits>

//#define GAUXC_DISABLE_TIMINGS

namespace GauXC {
namespace util  {

namespace detail {
  // TODO: write type_traits for GauXC
 
  template <class F, class... Args>
  struct has_void_return_type {

    // TODO: Write C++20 friendly version with removal of std::result_of
    static constexpr bool value = 
      std::is_same< std::invoke_result_t<F,Args...>, void >::value;

  };
  
}

class Timer {

  template <class Rep, class Period>
  using duration = std::chrono::duration<Rep,Period>;

  std::map< std::string, duration<double, std::milli>> timings_;

public:


  Timer()                               = default;
  Timer( const Timer& )                 = default;
  Timer( Timer&& ) noexcept             = default;
  Timer& operator=( const Timer& )      = default;
  Timer& operator=( Timer&& ) noexcept  = default;

  template <class OtherRep, class OtherPeriod>
  void add_timing( std::string name, duration<OtherRep,OtherPeriod> dur ) {
    timings_.insert_or_assign( name, duration<double,std::milli>(dur) );
  }

  template <class OtherRep, class OtherPeriod>
  void add_or_accumulate_timing( std::string name, 
                                 duration<OtherRep,OtherPeriod> dur ) {

    if( timings_.find( name ) != timings_.end() ) {
      timings_.at( name ) += dur;
    } else {
      add_timing( name, dur );
    }

  }

  template <typename Op>
  inline 
  std::enable_if_t< detail::has_void_return_type<Op>::value > 
    time_op( std::string name, const Op& op ) {

#ifndef GAUXC_DISABLE_TIMINGS
    auto st = std::chrono::high_resolution_clock::now();
    op();
    auto en = std::chrono::high_resolution_clock::now();

    duration< double, std::milli > dur( en - st );
    add_timing( name, dur );
#else
    op();
#endif

  }

  template <typename Op>
  inline 
  std::enable_if_t< not detail::has_void_return_type<Op>::value, 
                    std::invoke_result_t<Op>
                  > time_op( std::string name, const Op& op ) {
#ifndef GAUXC_DISABLE_TIMINGS
    auto st = std::chrono::high_resolution_clock::now();
    auto res = op();
    auto en = std::chrono::high_resolution_clock::now();

    duration< double, std::milli > dur( en - st );
    add_timing( name, dur );

    return res;
#else
    return op();
#endif
  }






  template <typename Op>
  inline 
  std::enable_if_t< detail::has_void_return_type<Op>::value > 
    time_op_accumulate( std::string name, const Op& op ) {

#ifndef GAUXC_DISABLE_TIMINGS
    auto st = std::chrono::high_resolution_clock::now();
    op();
    auto en = std::chrono::high_resolution_clock::now();

    duration< double, std::milli > dur( en - st );
    add_or_accumulate_timing( name, dur );
#else
    op();
#endif

  }

  template <typename Op>
  inline 
  std::enable_if_t< not detail::has_void_return_type<Op>::value, 
                    std::invoke_result_t<Op>
                  > time_op_accumulate( std::string name, const Op& op ) {

#ifndef GAUXC_DISABLE_TIMINGS
    auto st = std::chrono::high_resolution_clock::now();
    auto res = op();
    auto en = std::chrono::high_resolution_clock::now();

    duration< double, std::milli > dur( en - st );
    add_or_accumulate_timing( name, dur );

    return res;
#else
    return op();
#endif
  }





  template <class Rep = double, class Period = std::milli>
  inline duration<Rep,Period> get_duration( std::string name ) {
    return timings_.at(name);
  }

  inline const auto& all_timings() const { return timings_; }

};

}
}
