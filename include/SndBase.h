// SndBase.h 
// Base class and
// Binary operations
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

#ifndef _AURORA_SNDBASE_
#define _AURORA_SNDBASE_

#include <cstdint>
#include <numeric>
#include <vector>

namespace Aurora {
const int def_vsize = 64;
const double def_sr = 44100.;

/** SndBase class \n
    Aurora Library base class
*/
template <typename S> class SndBase {
protected:
  std::vector<S> sig;

public:
  /** Constructor \n
      vsize: signal vector size
  */
  SndBase(std::size_t vsize = def_vsize) : sig(vsize){};

  /** Vector size query \n
      returns the object vector size
  */
  std::size_t vsize() { return sig.size(); }

  /** Vector size setting \n
      n: new vector size
  */
  void vsize(std::size_t n) { sig.resize(n); }

  /** Vector access \n
      returns the object vector
  */
  const std::vector<S> &vector() { return sig; }

  /** Preallocate vector memory \n
      size: size of vector to reserve in memory \n
      this method does not change the vector size.
  */
  void prealloc(std::size_t size) { sig.reserve(size); }
};

/** BinOp class \n
    Binary operations
*/
template <typename S> class BinOp : public SndBase<S> {
  using SndBase<S>::sig;
  std::function<S(S, S)> op;

public:
  /** Constructor \n
      f: binary operation function \n
      vsize: signal vector size
  */
  BinOp(std::function<S(S, S)> f, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), op(f){};

  /** Binary operation \n
      a: scalar input \n
      s: signal input \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(S a, const std::vector<S> &s) {
    std::size_t n = 0;
    for (auto &out : sig)
      out = op(a, s[n++]);
    return sig;
  }

  /** Binary operation \n
      s: signal input \n
      a: scalar input \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &s, S a) {
    std::size_t n = 0;
    for (auto &out : sig)
      out = op(a, s[n++]);
    return sig;
  }

  /** Binary operation \n
      s1: signal input 1 \n
      s2: signal input 2 \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &s1,
                                   const std::vector<S> &s2) {
    std::size_t n = 0;
    for (auto &out : sig) {
      out = op(s1[n], s2[n]);
      n++;
    }
    return sig;
  }
};
} // namespace Aurora

#endif // _AURORA_SNDBASE_
