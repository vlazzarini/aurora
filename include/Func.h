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
    Generic templated function maps \n
    FN: function to be applied FN(arg, table)
    S: sample type
*/
template <typename S, S (*FN)(S, const std::vector<S> *)>
class Func : public SndBase<S> {
  using SndBase<S>::process;

public:
  /** Constructor \n
      vsize: signal vector size
  */
  Func(std::size_t vsize = def_vsize) : SndBase<S>(vsize){};

  /** Functional application \n
      in: input signal \n
      t: function table to be passed (optional) \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> *t = NULL) {
    std::size_t n = 0;
    auto fp = [&]() { return FN(in[n++], t); };
    return process(fp, in.size());
  }
};

} // namespace Aurora

#endif // _AURORA_FUNC_
