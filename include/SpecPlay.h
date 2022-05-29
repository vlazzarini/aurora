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

  /** SpecTable class \n
      spectral function table \n
      S - sample type
  */
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
  /** Constructor \n
       wn - analysis window \n
       hs - hopsize \n
       sr - sampling rate
  */
  SpecTable(const std::vector<S> &wn, std::size_t hsiz=def_hsize, S sr = def_sr)
    : win(wn), tab(0), hs(hsiz), fs(sr) {
    }

    /** Spectral table creation \n
        in - input data \n
        returns the number of spectral frames
     */
    std::size_t operator()(const std::vector<S> &in) {
      create(in);
      return tab.size();
    }

    /** Spectral table access \n
        returns the spectral table
    */
    const std::vector<std::vector<Aurora::specdata<S>>>
      &operator()() const {
      return tab;
    }

    /** returns the analysis hopsize */
    std::size_t hsize() const { return hs; }

    /** returns the number of frames in the table */
    std::size_t size() const { return tab.size(); }

    /** sets the sampling rate */
    void sr(S r) {
      fs = r;
    }

    /** clears the table memory */
    void clear() {
      tab.clear();
    }
  
  };

  /** SpecPlay class \n
      spectral function table player \n
      S - sample type
  */  
  template <typename S> class SpecPlay {
    SpecShift<S> shift;
    S sr;
    S rp;
    S shft, fscal;
    S bn, fine, tscal, beg, end, st;
    bool keep;


  public:
  /** Constructor \n
       sr - sampling rate \n
       fftsize - analysis window size
  */ 
  SpecPlay(S fs = def_sr, std::size_t fftsize = def_fftsize) : shift(fs,fftsize), sr(fs), rp(0),
      shft(0), fscal(1), bn(261), fine(1), tscal(1), beg(0), end(1), st(0), keep(0){ }

    /** play from start \n
        siz - spec table size
     */ 
    void onset(std::size_t siz) {
      rp = (st < end? st : end)*siz;
    }

    /** reset object \n 
       fs - sampling rate
    */
    void reset(S fs) {
      shift.reset(fs);
      sr  = fs;
      rp = 0;
    }

    /** set frequency shift (Hz) */
    void freqshift(S f) { shft = f; }

    /** set formant scaling factor */
    void formscal(S f) { fscal = f; }

    /** set base freq in Hz */
    void basefreq(S f) { bn = f; }

    /** set fine scaling factor */
    void finetune(S f) { fine = f; }

    /** set time scale factor */
    void timescale(S ts) { tscal = ts; }

    /** loop begin pos (0 - 1) */
    void loopbeg(S p) { beg = p >= 0 ? (p < 1 ? p : 1) : 0; }

    /** loop end pos (0 - 1) */
    void loopend(S p) { end = p >= 0 ? (p < 1 ? p : 1) : 0;; }

    /** start pos (0 - 1) */
    void start(S p) { st = p >= 0 ? (p < 1 ? p : 1) : 0;; }

    /** extract and keep formants */
    void keepform(bool b) { keep = b; }

    /** Spectral player \n
        samp - spectral table \n
        cps - frequency \n
        returns spectral frame \n
    */
    const std::vector<specdata<S>>
      &operator() (const SpecTable<S> &samp, S cps) {
      std::size_t siz = samp.size();
      if(siz == 0){
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
