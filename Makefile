CC = gcc 
flags = -Wall -pedantic -std=c99 -O2
libs = -lm -lfftw3 -lgmp -lflint
noflint = -lm -lfftw3

all: 
	cp src/sp_mwdc_flint_matcher.c.orig src/sp_mwdc_flint_matcher.c 
	cp src/sp_mwdc_flint_matcher.h.orig src/sp_mwdc_flint_matcher.h 
	make harness -B
	make library -B
	make examples -B
    
harness: src/sp_mwdc_harness.c src/sp_mwdc.c src/sp_mwdc_naive_matcher.c src/sp_mwdc_fft_matcher.c src/sp_mwdc_flint_matcher.c src/sp_mwdc_harness.h src/sp_mwdc.h src/sp_mwdc_naive_matcher.h src/sp_mwdc_fft_matcher.h src/sp_mwdc_flint_matcher.h src/sp_mwdc_fft_random_matcher.c src/sp_mwdc_fft_random_matcher.h
	cd src/  && $(CC) $(flags) sp_mwdc_harness.c sp_mwdc.c sp_mwdc_naive_matcher.c sp_mwdc_fft_matcher.c sp_mwdc_flint_matcher.c sp_mwdc_fft_random_matcher.c -o ../sp_mwdc_harness $(libs)
	cd src/  && $(CC) $(flags) sp_km_harness.c sp_km.c sp_km_naive_matcher.c sp_km_unbounded_matcher.c -o ../sp_km_harness $(libs)
	cd src/  && $(CC) $(flags) sp_create_input.c -o ../sp_create_input $(libs)


library: src/sp_mwdc_harness.c src/sp_mwdc.c src/sp_mwdc_naive_matcher.c src/sp_mwdc_fft_matcher.c src/sp_mwdc_flint_matcher.c src/sp_mwdc_harness.h src/sp_mwdc.h src/sp_mwdc_naive_matcher.h src/sp_mwdc_fft_matcher.h src/sp_mwdc_flint_matcher.h src/sp_mwdc_fft_random_matcher.c src/sp_mwdc_fft_random_matcher.h src/sp_km.c src/sp_km.h src/sp_km_naive_matcher.c src/sp_km_naive_matcher.h src/sp_km_unbounded_matcher.c src/sp_km_unbounded_matcher.h
	cd src/ && $(CC) $(flags) -c sp_mwdc.c
	cd src/ && $(CC) $(flags) -c sp_mwdc_naive_matcher.c
	cd src/ && $(CC) $(flags) -c sp_mwdc_flint_matcher.c
	cd src/ && $(CC) $(flags) -c sp_mwdc_fft_random_matcher.c sp_mwdc_fft_matcher.c
	cd src/ && $(CC) $(flags) -c sp_km.c
	cd src/ && $(CC) $(flags) -c sp_km_naive_matcher.c
	cd src/ && $(CC) $(flags) -c sp_km_unbounded_matcher.c
	ar rvs lib/libstringpedia.a src/*.o
	cp src/*.h include/

examples: examples/sp_mwdc_example.c examples/sp_km_example.c
	cd examples/ && $(CC) $(flags) sp_mwdc_example.c ../lib/libstringpedia.a -o sp_mwdc_example $(libs)
	cd examples/ && $(CC) $(flags) sp_km_example.c ../lib/libstringpedia.a -o sp_km_example $(libs)


noflint: 
	cp src/sp_mwdc_flint_matcher.c.blank src/sp_mwdc_flint_matcher.c 
	cp src/sp_mwdc_flint_matcher.h.blank src/sp_mwdc_flint_matcher.h
	make harness-noflint -B
	make library-noflint -B
	make examples-noflint -B

harness-noflint: src/sp_mwdc_harness.c src/sp_mwdc.c src/sp_mwdc_naive_matcher.c src/sp_mwdc_fft_matcher.c src/sp_mwdc_flint_matcher.c src/sp_mwdc_harness.h src/sp_mwdc.h src/sp_mwdc_naive_matcher.h src/sp_mwdc_fft_matcher.h src/sp_mwdc_flint_matcher.h src/sp_mwdc_fft_random_matcher.c src/sp_mwdc_fft_random_matcher.h
	cd src/  && $(CC) $(flags) sp_mwdc_harness.c sp_mwdc.c sp_mwdc_naive_matcher.c sp_mwdc_fft_matcher.c sp_mwdc_flint_matcher.c sp_mwdc_fft_random_matcher.c -o ../sp_mwdc_harness $(noflint)
	cd src/  && $(CC) $(flags) sp_km_harness.c sp_km.c sp_km_naive_matcher.c sp_km_unbounded_matcher.c -o ../sp_km_harness $(noflint)
	cd src/  && $(CC) $(flags) sp_create_input.c -o ../sp_create_input $(noflint)


library-noflint: src/sp_mwdc_harness.c src/sp_mwdc.c src/sp_mwdc_naive_matcher.c src/sp_mwdc_fft_matcher.c src/sp_mwdc_flint_matcher.c src/sp_mwdc_harness.h src/sp_mwdc.h src/sp_mwdc_naive_matcher.h src/sp_mwdc_fft_matcher.h src/sp_mwdc_flint_matcher.h src/sp_mwdc_fft_random_matcher.c src/sp_mwdc_fft_random_matcher.h src/sp_km.c src/sp_km.h src/sp_km_naive_matcher.c src/sp_km_naive_matcher.h src/sp_km_unbounded_matcher.c src/sp_km_unbounded_matcher.h
	cd src/ && $(CC) $(flags) -c sp_mwdc.c
	cd src/ && $(CC) $(flags) -c sp_mwdc_naive_matcher.c
	cd src/ && $(CC) $(flags) -c sp_mwdc_flint_matcher.c
	cd src/ && $(CC) $(flags) -c sp_mwdc_fft_random_matcher.c sp_mwdc_fft_matcher.c
	cd src/ && $(CC) $(flags) -c sp_km.c
	cd src/ && $(CC) $(flags) -c sp_km_naive_matcher.c
	cd src/ && $(CC) $(flags) -c sp_km_unbounded_matcher.c
	ar rvs lib/libstringpedia.a src/*.o
	cp src/*.h include/

examples-noflint: examples/sp_mwdc_example.c examples/sp_km_example.c
	cd examples/ && $(CC) $(flags) sp_mwdc_example.c ../lib/libstringpedia.a -o sp_mwdc_example $(noflint)
	cd examples/ && $(CC) $(flags) sp_km_example.c ../lib/libstringpedia.a -o sp_km_example $(noflint)



clean: clean-examples clean-objectfiles clean-lib clean-include

clean-examples: 
	rm -f sp_mwdc_harness
	rm -f sp_km_harness
	rm -f sp_create_input
	cd examples/ && rm -f sp_mwdc_example sp_km_example
	

clean-objectfiles: 
	cd src && rm -f *.o

clean-lib: 
	cd lib && rm -f *.a

clean-include:
	cd include && rm -f *.h
	

