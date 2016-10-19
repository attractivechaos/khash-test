#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "khash.h"
//KHASH_INIT(64, khint64_t, int, 1, kh_int64_hash_func2, kh_int64_hash_equal)
KHASH_MAP_INIT_INT64(64, int)
KHASH_MAP_INIT_STR(str, int)

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

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

char *read_seq(const char *fn, long *len)
{
	long l = 0, m = 0, m_buf = 0, l_buf;
	char *s = 0, *buf = 0;
	FILE *fp;
	*len = 0;
	fp = fn && strcmp(fn, "-")? fopen(fn, "r") : stdin;
	if (fp == 0) return 0;
	while ((l_buf = kgetline(&buf, &m_buf, fp)) != EOF) {
		if (l + l_buf >= m) {
			m = l + l_buf + 1;
			kroundup32(m);
			s = (char*)realloc(s, m);
		}
		strncpy(s + l, buf, l_buf);
		l += l_buf;
	}
	fclose(fp);
	s = (char*)realloc(s, l + 1);
	s[l] = 0;
	*len = l;
	return s;
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

int count_64(const char *seq, int len, int k)
{
	khash_t(64) *h;
    int i, l, size;
    uint64_t x, mask = (1ULL<<k*2) - 1;
	h = kh_init(64);
    for (i = l = 0, x = 0; i < len; ++i) {
        int absent, c = (uint8_t)seq[i] < 128? seq_nt4_table[(uint8_t)seq[i]] : 4;
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & mask;
            if (++l >= k) { // we find a k-mer
                khint_t itr;
                itr = kh_put(64, h, x, &absent); // only add one strand!
                if (absent) kh_val(h, itr) = 0;
                ++kh_val(h, itr);
            }
        } else l = 0, x = 0; // if there is an "N", restart
    }

	size = kh_size(h);
	kh_destroy(64, h);
	return size;
}

int main(int argc, char *argv[])
{
	long len;
	char *seq;
	clock_t t;
	int k = 31, cnt;

	if (argc == 1) {
		fprintf(stderr, "Usage: count <in.seq> [k=%d]\n", k);
		return 1;
	}
	if (argc >= 3) k = atoi(argv[2]);

	seq = read_seq(argv[1], &len);
	t = clock();
	cnt = count_64(seq, len, k);
	fprintf(stderr, "%d\t%.3f\n", cnt, (double)(clock() - t) / CLOCKS_PER_SEC);

	free(seq);
	return 0;
}
