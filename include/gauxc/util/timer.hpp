#pragma once
#include <chrono>
#include <map>
#include <string>

namespace GauXC {
namespace util  {

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

  template <typename Op>
  inline auto time_op( std::string name, const Op& op ) {
    auto st = std::chrono::high_resolution_clock::now();
    op();
    auto en = std::chrono::high_resolution_clock::now();

    duration< double, std::milli > dur( en - st );
    add_timing( name, dur );

    return dur;
  }

  template <class Rep = double, class Period = std::milli>
  inline duration<Rep,Period> get_duration( std::string name ) {
    return timings_.at(name);
  }

  inline const auto& all_timings() const { return timings_; }

};

}
}
