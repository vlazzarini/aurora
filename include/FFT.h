// fft.h
// Fast Fourier transform functions
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

#ifndef _AURORA_FFT_
#define _AURORA_FFT_

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Aurora {

/**
   constant indicating packed FFT format
 */
const bool packed = true;

/**
   constant indicating forward FFT direction
*/
static bool forward = true;

/**
   constant indicating inverse FFT direction
*/
static bool inverse = false;

static inline uint32_t np2(uint32_t n) {
  uint32_t v = 2;
  while (v < n)
    v <<= 1;
  return v;
}

/** FFT class  \n
    Radix-2 fast Fourier transform \n
    S: sample type
*/
template <typename S> class FFT {
  std::vector<std::complex<S>> c;
  bool pckd;
  std::size_t sz;
  bool norm;

  void reorder(std::complex<S> *s) {
    const std::size_t N = sz;
    for (std::size_t i = 0, j = 0, m = 0; i < N; i++) {
      if (j > i) {
        std::swap(s[i], s[j]);
      }
      m = N / 2;
      while (m >= 2 && j >= m) {
        j -= m;
        m /= 2;
      }
      j += m;
    }
  }

public:
  /** Constructor \n
      N: transform size \n
      packed: complex data format in packed form (true) or not \n
      nm: normalisation mode (forward or inverse transforms)
  */
  FFT(std::size_t N, bool packd = packed, bool nm = forward)
      : c(packd ? np2(N) / 2 : np2(N) / 2 + 1), pckd(packd), sz(np2(N) / 2),
        norm(nm){ };

  std::size_t size() { return 2 * c.size(); }

  /** In-place complex-to-complex FFT \n
      data: data to be transformed  \n
      dir: true for forward operation, false for inverse \n
  */
  void transform(std::complex<S> *s, bool dir) {
    const std::size_t N = sz;
    std::complex<S> wp, w, even, odd;
    reorder(s);
    for (std::size_t n = 1, i = 0; n < N; n *= 2) {
      S o = dir == forward ? -M_PI / n : M_PI / n;
      wp.real(std::cos(o)), wp.imag(std::sin(o));
      w = 1.;
      for (std::size_t m = 0; m < n; m++) {
        for (std::size_t k = m; k < N; k += n * 2) {
          i = k + n;
          even = s[k];
          odd = w * s[i];
          s[k] = even + odd;
          s[i] = even - odd;
        }
        w *= wp;
      }
    }
    dir = norm ? dir : !dir;
    if (dir == forward) {
      for (std::size_t n = 0; n < N; n++)
        s[n] /= N;
    }
  }

  /** Real-to-complex forward FFT \n
      r: input vector (real)
      returns pointer to complex data containing the transform
  */
  const std::complex<S> *transform(const std::vector<S> &r) {
    using namespace std::complex_literals;
    const std::size_t N = sz;
    std::complex<S> wp, w = 1., even, odd;
    S o, zro, nyq;
    S *s = reinterpret_cast<S *>(c.data());
    std::fill(c.begin(), c.end(), std::complex<S>(0, 0));
    std::copy(r.begin(), r.end(), s);
    transform(c.data(), forward);
    zro = c[0].real() + c[0].imag();
    nyq = c[0].real() - c[0].imag();
    c[0].real(zro), c[0].imag(nyq);
    o = -M_PI / N;
    wp.real(std::cos(o)), wp.imag(std::sin(o));
    w *= wp;
    for (std::size_t i = 1, j = 0; i < N / 2; i++) {
      j = N - i;
      even = S(.5) * (c[i] + conj(c[j]));
      odd = std::complex<S>(.5i) * (conj(c[j]) - c[i]);
      c[i] = even + w * odd;
      c[j] = conj(even - w * odd);
      w *= wp;
    }
    if (!pckd) {
      c[N].real(c[0].imag());
      c[0].imag(0.);
      c[N].imag(0.);
    }
    return c.data();
  }

  /** Complex-to-real inverse FFT \n
      sp - input vector (complex) \n
      returns pointer to real data containing the transform
  */
  const S *transform(const std::vector<std::complex<S>> &sp) {
    using namespace std::complex_literals;
    const std::size_t N = sz;
    std::fill(c.begin(), c.end(), std::complex<S>(0, 0));
    std::complex<S> wp, w = 1., even, odd;
    S o, zro, nyq;
    S *s = reinterpret_cast<S *>(c.data());
    std::copy(sp.begin(), sp.end(), c.begin());
    if (pckd)
      zro = c[0].real()*0.5, nyq = c[0].imag()*0.5;
    else {
      zro = c[0].real()*0.5;
      nyq = c[N].real()*0.5;
    }
    c[0].real(zro + nyq), c[0].imag(zro - nyq);
    o = M_PI / N;
    wp.real(std::cos(o)), wp.imag(std::sin(o));
    w *= wp;
    for (std::size_t i = 1, j = 0; i < N / 2; i++) {
      j = N - i;
      even = S(.5) * (c[i] + conj(c[j]));
      odd = std::complex<S>(.5i) * (c[i] - conj(c[j]));
      c[i] = even + w * odd;
      c[j] = conj(even - w * odd);
      w *= wp;
    }
    transform(c.data(), inverse);
    return s;
  }

  const std::vector<std::complex<S>> &operator() (const std::vector<S> &r) {
    transform(r);
    return c;
  }

  const S *operator() (const std::vector<std::complex<S>> &sp)
  {
    return transform(sp);
  }

  const std::vector<std::complex<S>> &vector() const { return c; }
  
  const S *data() const { return reinterpret_cast<const S *>(c.data()); }
  
};
} // namespace Aurora

#endif // _AURORA_FFT_
