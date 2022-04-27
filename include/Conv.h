// Conv.h
// Partitioned convolution
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

#include "FFT.h"
#include "SndBase.h"
#include <complex>
#include <vector>

namespace Aurora {
  const std::size_t def_psize = 2048;
  const bool ola = 1;
  const bool ols = 0;
  /** IR class \n
      Impulse response table
  */
  template <typename S> class IR {
    std::vector<std::vector<std::complex<S>>> parts;

    void create(const std::vector<S> &s, std::size_t psize) {
      std::vector<S> buffer(psize);
      std::size_t n = 0, k = 0, cnt = s.size();
      FFT<S> fft(2 * psize, !packed, inverse);

      while (cnt) {
	cnt = cnt > psize ? psize : cnt;
	std::copy(s.begin() + n, s.begin() + n + cnt, buffer.begin());
	fft.transform(buffer);
	auto v = fft.vector();
	std::copy(v.begin(), v.end(), parts[k++].begin());
	std::fill(buffer.begin(), buffer.end(), 0);
	n += cnt;
	cnt = s.size() - n;
      }
    }

  public:
    /** Constructor \n
	s: impulse response \n
	psize: partition size
    */
  IR(const std::vector<S> &s, std::size_t psize = def_psize)
    : parts(ceil((float)s.size() / psize),
	    std::move(std::vector<std::complex<S>>(psize + 1)))

      {
	create(s, psize);
      }

    /** IR spectrum access \n
	returns a vector of impulse response partitions
    */
    const std::vector<std::vector<std::complex<S>>> &spectrum() const {
      return parts;
    }

    /** IR partition size \n
	returns the size of each partition
    */
    const std::size_t psize() const { return parts[0].size() - 1; }

    /** IR length \n
	returns the number of partitions
    */
    const std::size_t nparts() const { return parts.size(); }

    /** Reset the IR table \n
	s: impulse response \n
	psize: partition size
    */
    void reset(const std::vector<S> &s, std::size_t psize = def_psize) {
      parts.clear();
      parts.resize(ceil((float)s.size() / psize));
      for (auto &part : parts)
	part = std::move(std::vector<std::complex<S>>(psize + 1));
      create(s, psize);
    }
  };

  /** Conv class \n
      Partitioned convolution
  */
  template <typename S> class Conv : public SndBase<S> {
    using SndBase<S>::process;
    using SndBase<S>::vector;
    const IR<S> *ir;
  protected:
    std::vector<std::vector<std::complex<S>>> del;
    std::vector<std::vector<std::complex<S>>> del2;
    std::vector<std::complex<S>> mix;
    std::vector<S> inbuf;
    std::vector<S> inbuf2;
    std::vector<S> olabuf;
    std::size_t p, sn, psize;
    FFT<S> fft;
    bool meth;

    void convol(const std::vector<std::vector<std::complex<S>>> &in1,
		const std::vector<std::vector<std::complex<S>>> &in2,
		std::size_t pp) {
      auto dl = in1.begin() + pp;
      auto &part = in2;
      std::fill(mix.begin(), mix.end(), 0);
      for (auto prt = part.rbegin(); prt != part.rend(); prt++, dl++) {
	if (dl == in1.end())
	  dl = in1.begin();
	auto dsamp = dl->begin();
	auto psamp = prt->begin();
	for (auto &mx : mix)
	  mx += (*dsamp++ * *psamp++);
      }
      fft.transform(mix);
    }

    S oladd(S in, S *bufin, const S *bufout, S *olabuf, std::size_t cnt,
	    std::size_t psz) {
      auto s = bufout[cnt] + olabuf[cnt];
      bufin[cnt] = in;
      olabuf[cnt] = bufout[cnt + psz];
      return s;
    }

    S olsave(S in, S *bufin, const S *bufout, S *unused, std::size_t p,
	     std::size_t psz) {
      (void)unused;
      auto s = bufout[p + psz];
      bufin[p] = bufin[p + psz];
      bufin[p + psz] = in;
      return s;
    }

    void transform(const std::vector<S> &in,
		   std::vector<std::vector<std::complex<S>>> &d) {
      fft.transform(in);
      auto &v = fft.vector();
      std::copy(v.begin(), v.end(), d[p].begin());
    }				   

  public:
    /** Constructor \n
	IR: impulse response table\n
	algo: convolution algorithm (ols = 0, ola = 1) \n
	vsize: vector size
    */
  Conv(const IR<S> *imp, bool algo = ols, std::size_t vsize = def_vsize)
    : SndBase<S>(vsize), ir(imp),
      del(imp->nparts(), std::vector<std::complex<S>>(imp->psize() + 1)),
      del2(0), mix(imp->psize() + 1), inbuf(2 * imp->psize()), inbuf2(0),
      olabuf(imp->psize()),
      p(0), sn(0), psize(imp->psize()), fft(imp->psize() * 2,
					    !packed, inverse), meth(algo){};

    /** Constructor \n
	len: convolution length \n
	psize: partition size \n
	vsize: vector size
    */
  Conv(std::size_t len,  std::size_t psiz = def_psize,
       std::size_t vsize = def_vsize)
    : SndBase<S>(vsize), ir(nullptr),
      del(std::ceil(len/psiz), std::vector<std::complex<S>>(psiz + 1)),
      del2(std::ceil(len/psiz), std::vector<std::complex<S>>(psiz + 1)),
      mix(psiz + 1), inbuf(2 * psiz), inbuf2(2 * psiz), olabuf(psiz),
      p(0), sn(0), psize(psiz), fft(psize * 2, !packed, inverse),
      meth(ola){};

    /** Convolution \n
	in: input \n
	scal: output ampltude scaling
    */
    const std::vector<S> &operator()(const std::vector<S> &in, S scal) {
      if(ir == nullptr) return vector();
      std::size_t n = 0;
      S *bufin = inbuf.data();
      S *obuff = olabuf.data();
      std::size_t sz = psize;
      return meth ? process(
			    [&]() {
			      auto s =
				oladd(in[n++], bufin, fft.data(),
				      obuff, sn, sz);
			      if (++sn == sz) {
				transform(inbuf,del);
				p = p == del.size() - 1 ? 0 : p + 1;			  
				convol(del,ir->spectrum(),p);
				std::fill(inbuf.begin()+psize,inbuf.end(),0);
				sn = 0;
			      }
			      return s * scal;
			    },
			    in.size())
	: process(
		  [&]() {
		    auto s =
		      olsave(in[n++], bufin, fft.data(),
			     obuff, sn, sz);
		    if (++sn == sz) {
		      transform(inbuf,del);
		      p = p == del.size() - 1 ? 0 : p + 1;			  
		      convol(del,ir->spectrum(),p);
		      sn = 0;
		    }
		    return s * scal;
		  },
		  in.size());
    }

    const std::vector<S> &operator()(const std::vector<S> &in1,
				     const std::vector<S> &in2, S scal) {
      if(del2.size() == 0) return vector();
      std::size_t n = 0;
      S *bufin = inbuf.data();
      S *bufin2 = inbuf2.data();
      S *obuff = olabuf.data();
      std::size_t sz = psize;
      return process( [&]() {
	  auto s =
	    oladd(in1[n], bufin, fft.data(), obuff, sn, sz);
	  bufin2[sn] = in2[n];
	  if (++sn == sz) {
	    transform(inbuf,del);
	    transform(inbuf2,del2);
	    p = p == del.size() - 1 ? 0 : p + 1;			  
	    convol(del,del2,p);
	    sn = 0;
	  }
	  n++;
	  return s *scal;
	},
	in1.size());
    }

    void reset(const IR<S> *imp) {
      ir = imp;
      psize = ir->psize();	
      del = std::vector<std::vector<std::complex<S>>>
	(ir->nparts(), std::vector<std::complex<S>>(psize + 1));
      mix = std::vector<std::complex<S>>(psize + 1);
      inbuf = std::vector<S>(2 * psize);
      olabuf = std::vector<S>(psize);
      p = 0;
      sn = 0;
      fft = FFT<S>(psize * 2, !packed, inverse);    
    }

    void reset(std::size_t len, std::size_t psiz = def_psize) {
      del = std::vector<std::vector<std::complex<S>>>
	(len/psize, std::vector<std::complex<S>>(psize + 1));
      del2 = std::vector<std::vector<std::complex<S>>>
	(len/psize,std::vector<std::complex<S>>(psize + 1));
      mix = std::vector<std::complex<S>>(psize + 1);
      inbuf = std::vector<S>(2 * psize);
      inbuf2 = std::vector<S>(2 * psize);
      olabuf = std::vector<S>(psize);
      p = 0;
      sn = 0;
      fft = FFT<S>(psize * 2, !packed, inverse);
      psize = psiz;     
      } 
  };

} // namespace Aurora
