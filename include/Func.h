// Func.h:
// Generic function maps
//
// (c) V Lazzarini, 2021
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE

#ifndef _AURORA_FUNC_
#define _AURORA_FUNC_

#include "SndBase.h"
#include <functional>

namespace Aurora {

/** Func class  \n
    Generic function maps \n
    S: sample type
*/
template <typename S> class Func : public SndBase<S> {
  using SndBase<S>::sig;
  std::function<S(S)> fun;

public:
  /** Constructor \n
      f: oscillator function \n
      fs: sampling rate \n
      vsize: signal vector size
  */
  Func(std::function<S(S)> f, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), fun(f){};

  /** Functional application \n
      in: input signal \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &in) {
    std::size_t n = 0;
    for (auto &s : sig)
      s = fun(in[n++]);
    return sig;
  }

  /** set the mapping function \n
     f: mapping function to be used
  */
  void func(std::function<S(S)> f) { fun = f; }
};
} // namespace Aurora

#endif // _AURORA_FUNC_
