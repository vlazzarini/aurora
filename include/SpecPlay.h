// SpecPlay.h
// Spectral Sample Player
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

#include "SpecStream.h"
#include "SpecShift.h"

namespace Aurora {

  template <typename S> 
    class SpecTable {
    const std::vector<S> &win;
    std::vector<std::vector<Aurora::specdata<S>>> tab;
    std::size_t N;
    std::size_t hs;
    S fs;

    void create(const std::vector<S> &in) {
      if(in.size() > 0) {
	Aurora::SpecStream<S> anal(win,hs,fs);
	tab.resize(in.size()/hs);
	std::size_t n = 0;
	std::vector<S> lfr(anal.hsize());
	for(auto &frame : tab) {
	  std::copy(in.data() + n,
		    in.data() + n + hs,
		    lfr.begin());
	  frame = anal(lfr);
	  n += hs;
	}
	if(n < in.size()) {
	  std::fill(lfr.begin(), lfr.end(), 0);
	  std::copy(in.begin() + n,
		    in.end(),
		    lfr.begin());
	  tab.push_back(anal(lfr));
	}
      }
      std::vector<float> amps(win.size()/2 + 1);
      for(auto &frame : tab) {
	std::size_t n = 0;
	for(auto &bin : frame) {
	  amps[n++] += bin.amp()/tab.size();
	}
      }
    }

  public:
  SpecTable(const std::vector<S> &wn, std::size_t hsiz=def_hsize, S sr = def_sr)
    : win(wn), tab(0), hs(hsiz), fs(sr) {
    }

    std::size_t operator()(const std::vector<S> &in) {
      create(in);
      return tab.size();
    }

    const std::vector<std::vector<Aurora::specdata<S>>>
      &operator()() const {
      return tab;
    }

    std::size_t hsize() const { return hs; }
    std::size_t fftsize() const { return win.size(); }
    std::size_t size() const { return tab.size(); }
    void sr(S r) {
      fs = r;
    }
    void clear() {
      tab.clear();
    }
  
  };

  
  template <typename S> class SpecPlay {
    SpecShift<S> shift;
    S sr;
    S rp;
    S siz, shft, fscal;
    S bn, fine, tscal, beg, end, st;
    bool keep;


  public:
  SpecPlay(S fs, std::size_t fftsize) : shift(fs,fftsize), sr(fs), rp(0), siz(0),
      shft(0), fscal(1), bn(261), fine(1), tscal(1), beg(0), end(1), st(0), keep(0){ }

    void onset() {
      rp = (st < end? st : end)*siz;
    }
    void reset(S fs) {
      shift.reset(fs);
      sr  = fs;
      rp = 0;
    }

    void size(std::size_t sz) { siz = sz; }
    std::size_t size() { return siz; }

    void freqshift(S f) { shft = f; }

    void formscal(S f) { fscal = f; }

    void basefreq(S f) { bn = f; }

    void finetune(S f) { fine = f; }

    void timescale(S ts) { tscal = ts; }

    void loopbeg(S t) { beg = t; }

    void loopend(S t) { end = t; }

    void start(S t) { st = t; }

    void keepform(bool b) { keep = b; }

    const std::vector<specdata<S>>
      &operator() (const SpecTable<S> &samp, S cps) {
      if(samp.size() != siz || samp.size() == 0){
	shift.reset(sr);
	return shift.frame();
      } 	
      shift.lock_formants(keep);
      shift(samp()[(int)rp],cps*fine/bn, shft, fscal);
      rp += tscal;
      if(end <= beg) beg = end;
      if(tscal >= 0) {
	rp = rp < end*siz ? rp : beg*siz;
      } else {
        while(rp < 0) rp += siz;
        rp = rp > beg*siz ? rp : end*siz - 1;
      }	 
      return shift.frame();
    }
  };

}
