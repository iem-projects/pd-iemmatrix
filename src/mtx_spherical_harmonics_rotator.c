/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *  this functions depends on the GNU scientific library
 *
 * Copyright (c) 2026, Thomas Röck
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */


#include <math.h>
#include <stdlib.h>
#include "iemmatrix.h"

#define SQRT2 1.4142135623730951

typedef struct _matrix_
{
    size_t size1;
    size_t size2;
    double *data;
} matrix;

matrix *matrix_calloc(size_t size1, size_t size2) {
    matrix *_matrix = (matrix *)malloc(sizeof(matrix));
    if (_matrix == NULL) {
        return NULL;
    }
    _matrix->size1 = size1;
    _matrix->size2 = size2;
    _matrix->data = (double *)calloc(size1 * size2, sizeof(double));
    if (_matrix->data == NULL) {
        free(_matrix);
        return NULL;
    }
    return _matrix;
}

double matrix_get(matrix *matrix, size_t i, size_t j)
{
    return matrix->data[i * matrix->size2 + j];
}

static void matrix_set_(matrix *matrix, size_t i, size_t j, double value)
{
    matrix->data[i * matrix->size2 + j] = value;
}

matrix *matrix_free_(matrix *matrix)
{
    if (matrix != NULL) {
        if (matrix->data != NULL) {
            free(matrix->data);
        }
        free(matrix);
    }
    return NULL;
}

typedef struct _ivanic_s_
{
    matrix *R;
    matrix *Rz_alpha;
    matrix *Rz_beta;
    matrix *Rz_gamma;
    matrix *ping;
    matrix *pong;
    size_t N;
} ivanic_s;


static void mat_mul_three(matrix *A, matrix *B, matrix *C, matrix *temp, matrix *result)
{
    size_t i, j, k;
    for (i = 0; i < A->size1; ++i)
    {
        for (j = 0; j < B->size2; ++j)
        {
            double sum = 0.0;
            for (k = 0; k < A->size2; ++k)
            {
                sum += matrix_get(A, i, k) * matrix_get(B, k, j);
            }
            matrix_set_(temp, i, j, sum);
        }
    }
    for (i = 0; i < temp->size1; ++i)
    {
        for (j = 0; j < C->size2; ++j)
        {
            double sum = 0.0;
            for (k = 0; k < temp->size2; ++k)
            {
                sum += matrix_get(temp, i, k) * matrix_get(C, k, j);
            }
            matrix_set_(result, i, j, sum);
        }
    }
}


static void fill_Rz(matrix *Rz, double cosx, double sinx)
{
    matrix_set_(Rz, 0, 0, cosx);
    matrix_set_(Rz, 0, 2, sinx);
    matrix_set_(Rz, 1, 1, 1.0);
    matrix_set_(Rz, 2, 0, -sinx);
    matrix_set_(Rz, 2, 2, cosx);
}

static void fill_Ry(matrix *Ry, double cosx, double sinx)
{
    matrix_set_(Ry, 0, 0, 1.0);
    matrix_set_(Ry, 1, 1, cosx);
    matrix_set_(Ry, 1, 2, -sinx);
    matrix_set_(Ry, 2, 1, sinx);
    matrix_set_(Ry, 2, 2, cosx);
}

static ivanic_s *ivanic_s_new(size_t N)
{
    ivanic_s *s = (ivanic_s *)malloc(sizeof(ivanic_s));
    if (s == NULL)
        return NULL;

    s->N = N;
    int l = (N + 1) * (N + 1);
    s->R = matrix_calloc(l, l);
    s->Rz_alpha = matrix_calloc(3, 3);
    s->Rz_beta = matrix_calloc(3, 3);
    s->Rz_gamma = matrix_calloc(3, 3);
    s->ping = matrix_calloc(3, 3);
    s->pong = matrix_calloc(3, 3);

    matrix *matrices[] = {s->R, s->Rz_alpha, s->Rz_beta, s->Rz_gamma, s->ping, s->pong};
    for (int i = 0; i < 6; ++i) {
        if (matrices[i] == NULL) {
            for (int j = 0; j < i; ++j)
                matrix_free_(matrices[j]);
        free(s);
        return NULL;
    }
    }

    return s;
}

static void ivanic_s_free(ivanic_s *s)
{
    if (s == NULL)
        return;

    matrix *matrices[] = {s->R, s->Rz_alpha, s->Rz_beta, s->Rz_gamma, s->ping, s->pong};
    for (int i = 0; i < 6; ++i) {
        if (matrices[i] != NULL)
            matrix_free_(matrices[i]);
        }
    free(s);
}

static double P(int i, int l, int mu, int m_prime, matrix *R, int R_lm1_offset)
{
    i += 2; // shift i by 2 to skip R0
    mu += l - 1 + R_lm1_offset;
    int m_prime_lm1 = m_prime + l - 1 + R_lm1_offset;
    int twolm2 = 2 * l - 2 + R_lm1_offset;
    if (abs(m_prime) < l)
    {
        return matrix_get(R, i, 2) * matrix_get(R, mu, m_prime_lm1);
    }
    if (m_prime == l)
    {
        return matrix_get(R, i, 3) * matrix_get(R, mu, twolm2) - matrix_get(R, i, 1 + 0) * matrix_get(R, mu, R_lm1_offset);
    }
    return matrix_get(R, i, 3) * matrix_get(R, mu, R_lm1_offset) + matrix_get(R, i, 1 + 0) * matrix_get(R, mu, twolm2);
}

static double U(int l, int m1, int m2, matrix *R, int R_lm1_offset)
{
    return P(0, l, m1, m2, R, R_lm1_offset);
}

static double V(int l, int m1, int m2, matrix *R, int R_lm1_offset)
{
    if (m1 == 0)
    {
        return P(1, l, 1, m2, R, R_lm1_offset) + P(-1, l, -1, m2, R, R_lm1_offset);
    }
    if (m1 == 1)
    {
        return SQRT2 * P(1, l, 0, m2, R, R_lm1_offset);
    }
    if (m1 == -1)
    {
        return SQRT2 * P(-1, l, 0, m2, R, R_lm1_offset);
    }
    if (m1 > 0)
    {
        return P(1, l, m1 - 1, m2, R, R_lm1_offset) - P(-1, l, -m1 + 1, m2, R, R_lm1_offset);
    }
    else
    {
        return P(1, l, m1 + 1, m2, R, R_lm1_offset) + P(-1, l, -m1 - 1, m2, R, R_lm1_offset);
    }
}

static double W(int l, int m1, int m2, matrix *R, int R_lm1_offset)
{
    if (m1 > 0)
    {
        return P(1, l, m1 + 1, m2, R, R_lm1_offset) + P(-1, l, -m1 - 1, m2, R, R_lm1_offset);
    }
    if (m1 < 0)
    {
        return P(1, l, m1 - 1, m2, R, R_lm1_offset) - P(-1, l, -m1 + 1, m2, R, R_lm1_offset);
    }
    return NAN;
}

static void uvw(double buffer[3], int l, int m1, int m2)
{
    static double denom;
    static const double one_half = 0.5;
    int abs_m1 = abs(m1);
    int m1_is_zero = (m1 == 0);
    if (abs(m2) < l)
    {
        denom = (double)((l + m2) * (l - m2));
    }
    else
    {
        denom = (double)((l << 1) * ((l << 1) - 1));
    }

    int u_num = (l + m1) * (l - m1);
    int v_num = (1 + m1_is_zero) * (l + abs_m1 - 1) * (l + abs_m1);
    int w_num = (l - abs_m1 - 1) * (l - abs_m1);

    buffer[0] = sqrt((double)u_num / denom);
    buffer[1] = one_half * sqrt((double)v_num / denom) * (double)(1 - 2 * m1_is_zero);
    buffer[2] = -one_half * sqrt((double)w_num / denom) * (double)(1 - m1_is_zero);
}

void ivanic_ruedenberg_rotation_matrix(ivanic_s *s, double alpha, double beta, double gamma)
{
    static double uvw_buffer[3];

    double cos_alpha = cos(alpha);
    double sin_alpha = sin(alpha);
    double cos_beta = cos(beta);
    double sin_beta = sin(beta);
    double cos_gamma = cos(gamma);
    double sin_gamma = sin(gamma);

    // Rz_alpha
    fill_Rz(s->Rz_alpha, cos_alpha, sin_alpha);

    // Rz_beta
    fill_Ry(s->Rz_beta, cos_beta, sin_beta);

    // Rz_gamma
    fill_Rz(s->Rz_gamma, cos_gamma, sin_gamma);

    // TODO R0
    // gsl_matrix_view R0_view = gsl_matrix_submatrix(s->R, 1, 1, 3, 3);
    // mat_mul_three(s->Rz_gamma, s->Ry_beta, s->Rz_alpha, s->temp, &(R0_view.matrix));
    // TODO copy into s

    int n_start, n_harmonics;

    // initialize matrices
    my_matrix_set(s->R, 0, 0, 1.0);

    if (s->N <= 1)
    {
        return;
    }

    int R_lm1_offset = 1;

    int m1, m2, nm1, nm2;
    double u_val, v_val, w_val;
    return;
    for (size_t n = 2; n <= s->N; ++n)
    {
        n_start = n * n + n;
        n_harmonics = (2 * n) + 1;
        R_lm1_offset = (n - 1) * (n - 1);
        for (m1 = -n; m1 <= n; ++m1)
        {
            for (m2 = -n; m2 <= n; ++m2)
            {
                u_val = 0.0;
                v_val = 0.0;
                w_val = 0.0;
                uvw(uvw_buffer, n, m1, m2);

                // ------ U ------
                if (fabs(uvw_buffer[0]) > 1e-15)
                {
                    // TODO compute_function_value(U, uvw_buffer[0], n, m1, m2, n_start, R, &M, P, mode);
                }
                // ------- V ------
                if (fabs(uvw_buffer[1]) > 1e-15)
                {
                    // TODO compute_function_value(V, uvw_buffer[1], n, m1, m2, n_start, R, &M, P, mode);
                }
                // ------- W ------
                if (fabs(uvw_buffer[2]) > 1e-15)
                {
                    // TODO compute_function_value(W, uvw_buffer[2], n, m1, m2, n_start, R, &M, P, mode);
                }
            }
        }
    }
}

typedef struct _mtx_spherical_harmonics_rotator_
{
    t_object x_obj;
    t_outlet *list_rot_out;
    t_atom *list_rot;
    ivanic_s *ivanic;
    double alpha;
    double beta;
    double gamma;
} t_mtx_spherical_harmonics_rotator;

static t_class *mtx_spherical_harmonics_rotator_class;

static t_mtx_spherical_harmonics_rotator *mtx_spherical_harmonics_rotator_new(t_symbol *s, int argc, t_atom *argv)
{
    int N;
    t_mtx_spherical_harmonics_rotator *x = (t_mtx_spherical_harmonics_rotator *)pd_new(mtx_spherical_harmonics_rotator_class);
    x->list_rot_out = outlet_new(&x->x_obj, gensym("matrix"));
    if (argc < 1)
    {
        N = 1;
    }
    if (argc >= 1)
    {
        N = (int)atom_getfloat(argv);
        if (N < 0)
        {
            N = 0;
        }
        x->ivanic = ivanic_s_new((size_t)N);
    }
    x->list_rot = (t_atom *)calloc((size_t)pow((N + 1), 4) + 2, sizeof(t_atom));
    x->alpha = 0.0;
    x->beta = 0.0;
    x->gamma = 0.0;
    return (void *)x;
}

static void mtx_spherical_harmonics_rotator_free(t_mtx_spherical_harmonics_rotator *x)
{
    if (x->list_rot)
    {
        free(x->list_rot);
    }
    if (x->ivanic)
    {
        ivanic_s_free(x->ivanic);
    }
}

static void mtx_spherical_harmonics_rotator_bang(t_mtx_spherical_harmonics_rotator *x)
{
    if (x->list_rot != 0)
    {
        outlet_anything(x->list_rot_out, gensym("matrix"),
                        (size_t)pow((x->ivanic->N + 1), 4) + 2, x->list_rot);
    }
}

static void mtx_spherical_harmonics_rotator_list(t_mtx_spherical_harmonics_rotator *x, t_symbol *s,
                                                  int argc, t_atom *argv)
{
    size_t num_harmonics, n;
    if (argc >= 1)
    {
        x->alpha = atom_getfloat(argv + 0);
    }
    if (argc >= 2)
    {
        x->beta = atom_getfloat(argv + 1);
    }
    if (argc >= 3)
    {
        x->gamma = atom_getfloat(argv + 2);
    }
    ivanic_ruedenberg_rotation_matrix(x->ivanic, x->alpha, x->beta, x->gamma);
    num_harmonics = (size_t)((x->ivanic->N + 1) * (x->ivanic->N + 1));
    SETFLOAT(x->list_rot, (t_float)num_harmonics);
    SETFLOAT(x->list_rot + 1, (t_float)num_harmonics);
    for (n = 0; n < num_harmonics * num_harmonics; ++n)
    {
        SETFLOAT(x->list_rot + n + 2, (t_float)x->ivanic->R->data[n]);
    }
    mtx_spherical_harmonics_rotator_bang(x);
}

void mtx_spherical_harmonics_rotator_setup(void)
{
    mtx_spherical_harmonics_rotator_class = class_new(
        gensym("mtx_spherical_harmonics_rotator"),
        (t_newmethod)mtx_spherical_harmonics_rotator_new,
        (t_method)mtx_spherical_harmonics_rotator_free,
        sizeof(t_mtx_spherical_harmonics_rotator),
        CLASS_DEFAULT,
        A_GIMME, 0);
    class_addbang(mtx_spherical_harmonics_rotator_class, (t_method)mtx_spherical_harmonics_rotator_bang);
    class_addlist(mtx_spherical_harmonics_rotator_class, (t_method)mtx_spherical_harmonics_rotator_list);
}


void iemtx_spherical_harmonics_rotator_setup(void)
{
    mtx_spherical_harmonics_rotator_setup();
}
