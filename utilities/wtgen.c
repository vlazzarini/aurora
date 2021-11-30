// wtgen.c
// Soundfile to C++ header wavetable
// depends on libsndfile
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

#include <sndfile.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, const char *argv[]) {
  if (argc > 6) {
    SNDFILE *fp;
    FILE *fpo;
    SF_INFO sfinfo;
    float f[BUFSIZ];
    int n, i, toread;
    size_t cnt, ends, starts, remain;
    fp = sf_open(argv[1], SFM_READ, &sfinfo);
    fpo = fopen(argv[2], "w");
    starts = atoi(argv[4]), ends = atoi(argv[5]);
    sf_seek(fp, starts, SEEK_SET);
    cnt = starts;
    remain = ends - starts;
    fprintf(fpo, "/* Wavetable created by %s */\n", argv[0]);
    fprintf(fpo, "/* source: %s (%lu - %lu samples) */\n", argv[1], starts,
            ends);
    fprintf(fpo, "#ifndef _%s_ \n", argv[3]);
    fprintf(fpo, "#define _%s_ \n", argv[3]);
    fprintf(fpo, "#include <vector>\n");
    fprintf(fpo, "namespace %s {\n", argv[3]);
    fprintf(fpo, "const std::vector<float> wave({\n");
    do {
      toread = remain > BUFSIZ ? BUFSIZ : remain;
      n = sf_readf_float(fp, f, toread);
      if (n == 0)
        break;
      cnt += n;
      remain -= n;
      for (i = 0; i < n * sfinfo.channels; i += sfinfo.channels) {
        fprintf(fpo, "%f,\n", f[i]);
      }
    } while (remain);
    fprintf(fpo, "});\n");
    fprintf(fpo, "/* base frequency */\n");
    fprintf(fpo, "const float base = %f;\n", atof(argv[6]));
    fprintf(fpo, "/* sampling rate */\n");
    fprintf(fpo, "const float fs = %f;\n", (float)sfinfo.samplerate);
    fprintf(fpo, "/* frequency ratio */\n");
    fprintf(fpo, "const float ratio = 1/(base*wave.size()/fs);\n");
    fprintf(fpo, "}\n");
    fprintf(fpo, "#endif /*_%s_ */\n", argv[3]);
    sf_close(fp);
    fclose(fpo);

    if (remain)
      printf("could not fully read requested range: end of file reached\n");
    printf("Wrote %lu samples (%f secs) to wavetable\n", cnt - starts,
           (float)(cnt - starts) / sfinfo.samplerate);
  } else
    printf("usage: %s infile outfile namespace start(samples) end(samples) "
           "basefr\n",
           argv[0]);

  return 0;
}
