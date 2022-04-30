// SpecBase.h
// Spectral Base class
//
// (c) V Lazzarini, 2022
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

#ifndef _AURORA_SPECBASE_
#define _AURORA_SPECBASE_

#include "SndBase.h"
#include "FFT.h"
#include <complex>
#include <cmath>


namespace Aurora {

  const std::size_t def_fftsize = 1024;
  const std::size_t def_hsize = 256;

  template <typename S>
    S unwrap(S ph) { return ph >= M_PI ? ph - twopi : (ph < -M_PI ? ph + twopi : ph); } 

  /** Spectral Data Type \n
      S: sample type
   */
  template <typename S> class specdata {
    std::complex<S>  bin;

  public:

  /** Spectral Data Type \n
      default constructor
   */
  specdata() : bin(0,0) { };

   /** Spectral Data Type
       amp - amplitude \n
       freq - freq
   */
  specdata(S amp, S freq) : bin(amp, freq) {};

   /** Spectral Data Type \n
       c - complex data in rectangular format
   */
  specdata(std::complex<S> c) : bin(c) {
      bin.real(std::abs(c));
      bin.imag(std::arg(c));
  }

  /** Cast operator \n
      returns complex data in rectangular format
  */ 
  operator std::complex<S> () const {
    return std::polar(bin.real(),bin.imag());
  }

  /** spectral data amplitude
   */
  S amp() const { return bin.real(); }

  /** spectral data frequency
   */
  S freq() const { return bin.imag(); }

  /** set spectral data amplitude
   */
  void amp(S a) { bin.real(a); }

  /** set spectral data frequency
   */
  void freq(S f) { bin.imag(f); }

  /** spectral data amplitude scaling
   */
  specdata &operator*=(S amp) {
      this->bin.real(amp*this->bin.real());
      return *this;
    }  
  
   /** spectral data amplitude scaling
   */
    specdata operator*(S amp) const {
      specdata res(*this);
      return res *= amp;
    }

   /** spectral data product
   */    
    friend specdata<S> operator*(S amp, const specdata<S> &d) {
      return d * amp;
    }

    /** phase difference \n
        oph - previous phase \n
        returns previous phase in object
    */
    S diff(S oph) {
      S ph = bin.imag();
      bin.imag(unwrap(ph-oph));
      return ph;
    }

     /** phase integration \n
        ph - current phase \n
        returns current phase in object
    */
     S integ(S ph) {
      bin.imag(ph+bin.imag());
      return bin.imag();
    }

     /** radians to Hz conversion \n
         cf - bin centre frequency in Hz \n
         fac - fs/(twopi*hopsize) conversion factor \n
         returns frequency in Hz.
     */
    S tocps(S cf, S fac) {
      return cf + bin.imag()*fac;
    };

     /** Hz to radian conversion \n
         cf - bin centre frequency in Hz \n
         fac - (twopi*hopsize)/fs conversion factor \n
         returns frequency in Hz.
     */    
    S fromcps(S cf, S fac) const {
      return (bin.imag() - cf)*fac;
    }
  
  };

  /** SpecBase class \n
     Spectral base class \n
     S: sample type
  */  
  template <typename S>
   class SpecBase { 
      std::vector<specdata<S>> spec;
      std::size_t hs;
      std::size_t fcnt;

    protected:
      std::vector<specdata<S>> &get_spec(){ return spec; }
      void fcount_incr() { fcnt++; }

    public:
        /** Constructor \n
        size: spectral frame size
        hsize: stream hopsize
        */
      SpecBase(std::size_t size = def_fftsize) :
      spec(size/2 + 1), fcnt(0) { };

      /** spectral frame size 
       */
      std::size_t size() const { return (spec.size() - 1)*2; }


      /** stream framecount
       */
      std::size_t framecount() const { return fcnt; }

      /** Frame access \n
         returns the object frame
      */
      const std::vector<specdata<S>> &frame() const { return spec; }
  
  };

    /** Ceps class \n
      Spectral envelope extraction
  */
  template <typename S>
    class Ceps : public SndBase<S>  {
    using SndBase<S>::get_sig;
    std::vector<std::complex<S>> spec;
    FFT<S> fft;   

   public:
   /** Constructor \n
       sixe: spectral frame size 
    */   
   Ceps(std::size_t size = def_fftsize) : SndBase<S>(size/2 + 1), spec(size/4 + 1),
      fft(size/2, !packed) { };

    /** Spectral envelope extraction \n
        in: input spectral frame \n
        coefs: number of cepstral coeffs retained \n
        returns the spectral envelope (non-negative frequencies only)
     */
    const std::vector<S> &operator()(const std::vector<specdata<S>> &in, std::size_t coefs){
      std::size_t n = 0;
      auto &mags = get_sig();
      std::transform(in.begin(), in.end(), mags.begin(),
		     [](const specdata<S> &in) { return in.amp() >= 0 ? std::log(in.amp()) : 0;} );  
      auto s = fft.transform(mags);
      std::fill(spec.begin(),spec.end(),0);
      std::copy(s,s+coefs,spec.begin());
      auto w = fft.transform(spec);
      for(auto &m : mags) m = std::exp(w[n++]);
      return mags;
    }
  };  
}

#endif
