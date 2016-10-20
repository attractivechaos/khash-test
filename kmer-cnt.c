#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>

#include "khash.h"
//KHASH_INIT(64, khint64_t, int, 1, kh_int64_hash_func2, kh_int64_hash_equal)
KHASH_MAP_INIT_INT64(64, int)
KHASH_MAP_INIT_STR(str, int)

typedef struct {
	char *s;
	khint_t h;
} cache_t;

#define cache_eq(a, b) ((a).h == (b).h && strcmp((a).s, (b).s) == 0)
#define cache_hash(a) ((a).h)
KHASH_INIT(cache, cache_t, int, 1, cache_hash, cache_eq)

long kgetline(char **buf, long *m_buf, FILE *fp)
{
    long tot = 0, max = 0;
    char *p;
    if (*m_buf == 0) { // empty buffer; allocate
        *m_buf = 256;   // initial size; could be larger
        *buf = (char*)malloc(*m_buf); // FIXME: check NULL
    }
    for (p = *buf, max = *m_buf;;) {
        long l, old_m;
        if (fgets(p, max, fp) == NULL)
            return tot? tot : EOF; // reach end-of-file
        for (l = 0; l < max; ++l)
            if (p[l] == '\n') {
                p[l] = 0;
                return tot + l;
            }
        old_m = *m_buf;
        *m_buf <<= 1; // incomplete line; double the buffer
        *buf = (char*)realloc(*buf, *m_buf); // check NULL
        max = (*m_buf) - old_m;
        p = (*buf) + old_m - 1; // point to the end of partial line
    }
    return tot;
}

unsigned char seq_nt4_table[128] = { // Table to change "ACGTN" to 01234
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int count_64(const char *fn, int k, int min_occ)
{
	khash_t(64) *h;
    long l_buf, m_buf = 0, size;
	khint_t itr;
	char *buf = 0;
	FILE *fp;

	if ((fp = fn && strcmp(fn, "-")? fopen(fn, "r") : stdin) == 0) return -1;
	h = kh_init(64);
	while ((l_buf = kgetline(&buf, &m_buf, fp)) != EOF) {
    	uint64_t x, mask = (1ULL<<k*2) - 1;
		int i, l;
		for (i = l = 0, x = 0; i < l_buf; ++i) {
			int absent, c = (uint8_t)buf[i] < 128? seq_nt4_table[(uint8_t)buf[i]] : 4;
			if (c < 4) { // not an "N" base
				x = (x << 2 | c) & mask;
				if (++l >= k) { // we find a k-mer
					khint_t itr;
					itr = kh_put(64, h, x, &absent);
					if (absent) kh_val(h, itr) = 0;
					++kh_val(h, itr);
				}
			} else l = 0, x = 0; // if there is an "N", restart
		}
	}
	fclose(fp);
	free(buf);

	for (itr = 0, size = 0; itr != kh_end(h); ++itr)
		if (kh_exist(h, itr) && kh_val(h, itr) >= min_occ)
			++size;
			
	kh_destroy(64, h);
	return size;
}

int count_str(const char *fn, int k, int min_occ)
{
	khash_t(str) *h;
    long l_buf, m_buf = 0, size;
	khint_t itr;
	char *buf = 0, buf2[64];
	FILE *fp;

	if ((fp = fn && strcmp(fn, "-")? fopen(fn, "r") : stdin) == 0) return -1;
	buf2[k] = 0;
	h = kh_init(str);
	while ((l_buf = kgetline(&buf, &m_buf, fp)) != EOF) {
		int i, l;
		for (i = l = 0; i < l_buf; ++i) {
			int absent, c = (uint8_t)buf[i] < 128? seq_nt4_table[(uint8_t)buf[i]] : 4;
			if (c < 4) { // not an "N" base
				if (++l >= k) { // we find a k-mer
					khint_t itr;
					strncpy(buf2, &buf[i + 1 - k], k);
					itr = kh_put(str, h, buf2, &absent);
					if (absent) {
						kh_key(h, itr) = strdup(buf2);
						kh_val(h, itr) = 0;
					}
					++kh_val(h, itr);
				}
			} else l = 0; // if there is an "N", restart
		}
	}
	fclose(fp);
	free(buf);

	for (itr = 0, size = 0; itr != kh_end(h); ++itr)
		if (kh_exist(h, itr)) {
			if (kh_val(h, itr) >= min_occ) ++size;
			free((char*)kh_key(h, itr));
		}
			
	kh_destroy(str, h);
	return size;
}

int count_cache(const char *fn, int k, int min_occ)
{
	khash_t(cache) *h;
    long l_buf, m_buf = 0, size;
	khint_t itr;
	cache_t cc;
	char *buf = 0, buf2[64];
	FILE *fp;

	if ((fp = fn && strcmp(fn, "-")? fopen(fn, "r") : stdin) == 0) return -1;
	buf2[k] = 0;
	cc.s = buf2;
	h = kh_init(cache);
	while ((l_buf = kgetline(&buf, &m_buf, fp)) != EOF) {
		int i, l;
		for (i = l = 0; i < l_buf; ++i) {
			int absent, c = (uint8_t)buf[i] < 128? seq_nt4_table[(uint8_t)buf[i]] : 4;
			if (c < 4) { // not an "N" base
				if (++l >= k) { // we find a k-mer
					khint_t itr;
					strncpy(cc.s, &buf[i + 1 - k], k);
					cc.h = __ac_X31_hash_string(cc.s);
					itr = kh_put(cache, h, cc, &absent);
					if (absent) {
						kh_key(h, itr).h = cc.h;
						kh_key(h, itr).s = strdup(cc.s);
						kh_val(h, itr) = 0;
					}
					++kh_val(h, itr);
				}
			} else l = 0; // if there is an "N", restart
		}
	}
	fclose(fp);
	free(buf);

	for (itr = 0, size = 0; itr != kh_end(h); ++itr)
		if (kh_exist(h, itr)) {
			if (kh_val(h, itr) >= min_occ) ++size;
			free(kh_key(h, itr).s);
		}
			
	kh_destroy(cache, h);
	return size;
}

int main(int argc, char *argv[])
{
	clock_t t;
	int c, k = 31, min_occ = 3, size, is_6 = 0, is_s = 0, is_c = 0;

	while ((c = getopt(argc, argv, "k:6sc")) >= 0)
		if (c == 'k') k = atoi(optarg);
		else if (c == '6') is_6 = 1;
		else if (c == 's') is_s = 1;
		else if (c == 'c') is_c = 1;

	if (argc == optind) {
		fprintf(stderr, "Usage: kmer-cnt [-6s] [-k %d] <in.seq>\n", k);
		return 1;
	}

	if (is_6) {
		t = clock();
		size = count_64(argv[optind], k, min_occ);
		fprintf(stderr, "%d\t%.3f\n", size, (double)(clock() - t) / CLOCKS_PER_SEC);
	}
	if (is_s) {
		t = clock();
		size = count_str(argv[optind], k, min_occ);
		fprintf(stderr, "%d\t%.3f\n", size, (double)(clock() - t) / CLOCKS_PER_SEC);
	}
	if (is_c) {
		t = clock();
		size = count_cache(argv[optind], k, min_occ);
		fprintf(stderr, "%d\t%.3f\n", size, (double)(clock() - t) / CLOCKS_PER_SEC);
	}
	return 0;
}
