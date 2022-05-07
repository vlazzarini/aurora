// SpecStream.h
// Streaming spectral analysis and synthesis
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

#ifndef _AURORA_SPECSTREAM_
#define _AURORA_SPECSTREAM_

#include <algorithm>
#include "SpecBase.h"

namespace Aurora {


  /** SpecStream class \n
      Spectral stream \n
      S: sample type
  */    
  template <typename S>
    class SpecStream : public SpecBase<S> {
    using SpecBase<S>::get_spec;
    using SpecBase<S>::fcount_incr;
    std::size_t hs;
    std::vector<S> buf;
    std::vector<S> wbuf;
    std::vector<S> oph;
    const std::vector<S> &win;
    FFT<S> fft;
    S fac, c;
    std::size_t dm, hnum, pos;
 
    void analysis(const std::vector<S> &in) {
      std::size_t n = 0;
      auto &v = fft(in);
      for(auto &s : get_spec()) {
	s = v[n];
	oph[n] = s.diff(oph[n]);
	s.freq(s.tocps(n*c, fac));
	n++;
      }
    }
    
  public:
   using SpecBase<S>::size;
    /** Constructor \n
        window: short-time Fourier transform window \n
        hsize: stream hopsize \n
        fs: sampling rate
    */
  SpecStream(const std::vector<S> &window, std::size_t hsize = def_hsize,
	     S fs = def_sr) :
   SpecBase<S>(window.size()), hs(hsize), buf(window.size()), wbuf(window.size()),
      oph(window.size()/2 + 1), win(window), fft(window.size(), !packed),
      fac(fs/(twopi*hsize)), c(fs/window.size()), dm(window.size()/hsize),
      hnum(dm-1), pos(0) { };

    /** Spectral stream analysis \n
        in: signal input \n
        returns the current spectral frame (non=negative frequencies only)
    */
    const std::vector<specdata<S>> &operator()(const std::vector<S> &in) {
      std::size_t vsize = in.size();
      if(vsize > hs) vsize = hs;
      std::size_t samps = vsize + pos;
      if(samps > hs) samps = vsize - samps + hs; 
      else samps = vsize;
      std::copy(in.begin(), in.begin() + samps, buf.begin() + pos + hs*hnum);
      pos += samps;
      if(pos == hs) {
        std::size_t n = 0, offs = hs*(dm - hnum - 1);
	std::size_t N = win.size();
        for (auto &s : wbuf) {
	  s = buf[n]*win[(n+offs)%N];
	  n++;
	}  
	analysis(wbuf);
	pos = 0;
        hnum = hnum != dm - 1 ? hnum + 1 : 0;
        if(samps != vsize) { std::copy(in.begin() + samps, in.begin() + vsize,
				       buf.begin() + hs*hnum);
          pos += vsize - samps;
	}
	fcount_incr();
	}
      return get_spec();
    }

    
   /** spectral hopsize 
     */
   std::size_t hsize() const { return hs; }

    /** reset the stream parameters \n
        fs - sampling rate 
    */ 
    void reset(S fs) {
      auto &ss = get_spec();
      fac = fs/(twopi*hsize());
      c = fs/win.size();
      std::fill(buf.begin(),buf.end(),0);
      std::fill(ss.begin(),ss.end(),specdata<S>(0,0));
      pos = 0;
    }
      
  };

  /** SpecSynth class \n
      STFT resynthesis \n
      S: sample type
  */    
  template <typename S>
    class SpecSynth : public SndBase<S> {
    using SndBase<S>::get_sig;
    std::vector<std::vector<S>> buffers;
    std::vector<std::complex<S>> spec;
    std::vector<double> ph;
    const std::vector<S> &win;
    FFT<S> fft;
    std::size_t dm, hsize;
    std::vector<std::size_t> count;
    S fac, c;

    const S *synthesis(const std::vector<specdata<S>> &in) {
      specdata<S> bin;
      std::size_t n = 0;
      for(auto &s : spec) {
	bin = in[n];
        bin.freq(bin.fromcps(n*c, fac));
	ph[n] = bin.integ(ph[n]);
	s = bin;
	n++;
      }
      return fft(spec);
    }
    
  public:

    /** Constructor \n
	window: short-time Fourier transform window \n
	hsize: stream hopsize \n
	fs: sampling rate \n
	vsize: output vector size
    */  
  SpecSynth(const std::vector<S> &window, std::size_t hsiz = def_hsize,
	    S fs = def_sr, std::size_t vsize = def_vsize) :
    SndBase<S>(vsize), buffers(window.size()/hsiz,
			       std::vector<S>(window.size())),
      spec(window.size()/2 + 1), ph(window.size()/2 + 1), win(window),
      fft(window.size(), !packed), dm(window.size()/hsiz),
      hsize(hsiz), count(dm), fac(twopi*hsiz/fs),
      c(fs/window.size()) {
      std::size_t n = 1;      
      for (auto &cnt : count) cnt = (dm - n++)*hsize;
    };


    /** Spectral stream synthesis \n
        in: input spectral frame  \n
        returns the output signal vector
    */    
    const std::vector<S> &operator() (const std::vector<specdata<S>> &in) {  
      std::size_t size = win.size();
      for(auto &ss : get_sig()) {
	ss = 0;
	std::size_t j = 0;
	for (auto &cnt : count) { 
	  auto &buf = buffers[j];
	  ss += buf[cnt++];
	  if(cnt == size) {
	    std::size_t n = 0;
	    auto s = synthesis(in);
	    std::size_t offs = hsize*j;
	    for(auto &b : buf) {
	      b = s[(n+offs)%size]*win[n];
	      n++;
	    }
	    cnt = 0;
	  }
	  j++;
	}
      }
      return get_sig();
    }
    
    /** reset the stream parameters \n
        fs - sampling rate 
    */ 
    void reset(S fs) {
      auto &ss = get_sig();
      fac = twopi*hsize/fs;
      c = fs/win.size();
      std::size_t n = 1;      
      for(auto &buf : buffers) {
	std::fill(buf.begin(), buf.end(), 0);
	count[n-1] = (dm - n)*hsize;
	n++;
      }
      std::fill(spec.begin(), spec.end(), std::complex<S>(0,0));
      std::fill(ss.begin(), ss.end(), 0);
      std::fill(ph.begin(), ph.end(), 0);
    }
    
  };

}
#endif
