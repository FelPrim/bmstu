#include <string.h>
#include <assert.h>

static double det_small(const Matrix *m){
    unsigned short n = m->rows;
    assert(m->rows == m->columns);
    switch (n){
    	case 0:
	       	return 0.0;
	case 1:
	       	return m->data[0];
	case 2:
        	return m->data[0]*m->data[3]
		     - m->data[1]*m->data[2];
    	
	case 3:
    		return m->data[0]*m->data[4]*m->data[8]
    	      		+ m->data[1]*m->data[5]*m->data[6]
    	      		+ m->data[2]*m->data[3]*m->data[7]
    	      		- m->data[2]*m->data[4]*m->data[6]
    	      		- m->data[1]*m->data[3]*m->data[8]
    	      		- m->data[0]*m->data[5]*m->data[7];
    }
}

static int contains_u16(const unsigned short *arr, unsigned short cnt, unsigned short v){
    for (unsigned short i = 0; i < cnt; ++i) if (arr[i] == v) return 1;
    return 0;
}

/* find the smallest row index in [0..n-1] not present in used_rows[0..used_count-1] */
static unsigned short next_free_row(const unsigned short *used_rows, unsigned short used_count, unsigned short n){
    for (unsigned short i = 0; i < n; ++i) if (!contains_u16(used_rows, used_count, i)) return i;
    /* should not happen if used_count < n */
    return n; /* sentinel */
}

/* collect remaining (not-used) rows and cols into out arrays; out_count = n - used_count */
static void collect_remaining(const Matrix *m, const unsigned short *used_rows, const unsigned short *used_cols,
                              unsigned short used_count,
                              unsigned short *out_rows, unsigned short *out_cols, unsigned short *out_count)
{
    unsigned short n = m->rows;
    unsigned short rr = 0, rc = 0;
    for (unsigned short i = 0; i < n; ++i){
        if (!contains_u16(used_rows, used_count, i)) out_rows[rr++] = i;
    }
    for (unsigned short j = 0; j < n; ++j){
        if (!contains_u16(used_cols, used_count, j)) out_cols[rc++] = j;
    }
    assert(rr == rc);
    *out_count = rr;
}

/* Frame for iterative DFS */
typedef struct {
    unsigned short choose_row; /* row expanded at this level */
    unsigned short next_col;   /* next column index to try (0..n) */
    int selected_col;          /* currently selected column for this level, -1 if none */
    double partial_mult;       /* accumulated multiplier up to this frame (product of a_ij and signs above) */
} Frame;

/* Public function - iterative Laplace expansion */
double determinant(const Matrix *m){
    assert(m && m->rows == m->columns);
    unsigned short n = m->rows;
    if (n <= 3) return det_small(m);

    /* VLA buffers on stack (no malloc) */
    Frame frames[n];
    unsigned short used_rows[n];
    unsigned short used_cols[n];

    /* init */
    for (unsigned short i = 0; i < n; ++i) { used_rows[i] = used_cols[i] = (unsigned short)(-1); }
    unsigned short depth = 1;            /* number of frames on stack (frames[0] is root) */
    unsigned short used_count = 0;       /* number of selected (row,col) pairs currently active */

    frames[0].choose_row = 0;            /* first free row is 0 */
    frames[0].next_col = 0;
    frames[0].selected_col = -1;
    frames[0].partial_mult = 1.0;

    const double *A = m->data;
    double total = 0.0;

    while (depth > 0){
        Frame *f = &frames[depth-1];
        /* find next free column for this frame starting from f->next_col */
        unsigned short j;
        for (j = f->next_col; j < n; ++j) if (!contains_u16(used_cols, used_count, j)) break;

        if (j >= n){
            /* no more columns to try at this frame -> pop */
            depth--; /* pop */
            if (depth == 0) break; /* finished */
            /* when returning to parent, unmark last used pair that belonged to parent selection */
            assert(used_count > 0);
            used_count--;
            frames[depth-1].selected_col = -1; /* parent selection ended */
            continue;
        }

        /* mark that next time we should continue after this j */
        f->next_col = j + 1;

        /* select column j for this row */
        used_rows[used_count] = f->choose_row;
        used_cols[used_count] = j;
        used_count++;
        f->selected_col = (int)j;

        /* compute new multiplier */
        int sign = (((unsigned int)f->choose_row + (unsigned int)j) & 1) ? -1 : 1;
        double aij = A[f->choose_row * n + j];
        double new_mult = f->partial_mult * (double)sign * aij;

        unsigned short rem = n - used_count;
        if (rem == 0){
            /* 0x0 minor => determinant 1 */
            total += new_mult * 1.0;
            /* unmark this selection and continue iterating this frame */
            used_count--;
            f->selected_col = -1;
            continue;
        }

        if (rem <= 3){
            unsigned short rows_rem[3], cols_rem[3];
            unsigned short rcnt;
            collect_remaining(m, used_rows, used_cols, used_count, rows_rem, cols_rem, &rcnt);
            double sub;
            if (rcnt == 1){
                sub = A[ rows_rem[0]*n + cols_rem[0] ];
            } else if (rcnt == 2){
                double a = A[ rows_rem[0]*n + cols_rem[0] ];
                double b = A[ rows_rem[0]*n + cols_rem[1] ];
                double c = A[ rows_rem[1]*n + cols_rem[0] ];
                double d = A[ rows_rem[1]*n + cols_rem[1] ];
                sub = a*d - b*c;
            } else { /* rcnt == 3 */
                double a00 = A[ rows_rem[0]*n + cols_rem[0] ];
                double a01 = A[ rows_rem[0]*n + cols_rem[1] ];
                double a02 = A[ rows_rem[0]*n + cols_rem[2] ];
                double a10 = A[ rows_rem[1]*n + cols_rem[0] ];
                double a11 = A[ rows_rem[1]*n + cols_rem[1] ];
                double a12 = A[ rows_rem[1]*n + cols_rem[2] ];
                double a20 = A[ rows_rem[2]*n + cols_rem[0] ];
                double a21 = A[ rows_rem[2]*n + cols_rem[1] ];
                double a22 = A[ rows_rem[2]*n + cols_rem[2] ];
                sub = a00*a11*a22 + a01*a12*a20 + a02*a10*a21
                    - a02*a11*a20 - a01*a10*a22 - a00*a12*a21;
            }
            total += new_mult * sub;
            /* unmark selection and continue */
            used_count--;
            f->selected_col = -1;
            continue;
        }

        /* rem > 3 : push child frame */
        unsigned short child_row = next_free_row(used_rows, used_count, n);
        assert(child_row < n);
        frames[depth].choose_row = child_row;
        frames[depth].next_col = 0;
        frames[depth].selected_col = -1;
        frames[depth].partial_mult = new_mult;
        depth++;
        /* continue loop to process new top frame */
    }

    return total;
}
