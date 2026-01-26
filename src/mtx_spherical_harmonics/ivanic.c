#include "ivanic.h"
#include "math.h"
#include "iemmatrix.h"
#include "iemmatrix_stub.h"
#ifdef HAVE_LIBGSL
#include <gsl/gsl_linalg.h>
#else
#include "stub/gsl.h"
#endif

IEMMATRIX_DECLARE_ALLOCFREE2_STUB(my_matrix);
#define SQRT2 1.4142135623730951
int have_gsl = 0;

typedef struct _ivanic_s_
{
    gsl_matrix *R; // Rotation matrix
    gsl_matrix *Rz_alpha;
    gsl_matrix *Ry_beta;
    gsl_matrix *Rz_gamma;
    gsl_matrix *temp;
    size_t N; // Maximum degree
} ivanic_s;

// static void mat_mul_three(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C, gsl_matrix *temp, gsl_matrix *result)
// {
//     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, temp);
//     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, C, 0.0, result);
// }

static void fill_Rz(gsl_matrix *Rz, double cosx, double sinx)
{
    my_matrix_set(Rz, 0, 0, cosx);
    my_matrix_set(Rz, 0, 2, sinx);
    my_matrix_set(Rz, 1, 1, 1.0);
    my_matrix_set(Rz, 2, 0, -sinx);
    my_matrix_set(Rz, 2, 2, cosx);
}

static void fill_Ry(gsl_matrix *Ry, double cosx, double sinx)
{
    my_matrix_set(Ry, 0, 0, 1.0);
    my_matrix_set(Ry, 1, 1, cosx);
    my_matrix_set(Ry, 1, 2, -sinx);
    my_matrix_set(Ry, 2, 1, sinx);
    my_matrix_set(Ry, 2, 2, cosx);
}

ivanic_s *ivanic_s_new(size_t N)
{
    post("ivanic_s_new called with N=%d", (int)N);
    post("have_gsl=%d", have_gsl);
    ivanic_s *s = (ivanic_s *)malloc(sizeof(ivanic_s));

    if (s == NULL)
    {
        return NULL;
    }
    s->N = N;
    int l = (N + 1) * (N + 1);
    s->R = (gsl_matrix *)my_matrix_calloc(l, l);
    if (s->R == NULL)
    {
        free(s);
        return NULL;
    }
    s->Rz_alpha = my_matrix_calloc(3, 3);
    if (s->Rz_alpha == NULL)
    {
        my_matrix_free(s->R);
        free(s);
        return NULL;
    }
    s->Ry_beta = my_matrix_calloc(3, 3);
    if (s->Ry_beta == NULL)
    {
        my_matrix_free(s->Rz_alpha);
        my_matrix_free(s->R);
        free(s);
        return NULL;
    }
    s->Rz_gamma = my_matrix_calloc(3, 3);
    if (s->Rz_gamma == NULL)
    {
        my_matrix_free(s->Ry_beta);
        my_matrix_free(s->Rz_alpha);
        my_matrix_free(s->R);
        free(s);
        return NULL;
    }
    s->temp = my_matrix_calloc(3, 3);
    if (s->temp == NULL)
    {
        my_matrix_free(s->Rz_gamma);
        my_matrix_free(s->Ry_beta);
        my_matrix_free(s->Rz_alpha);
        my_matrix_free(s->R);
        free(s);
        return NULL;
    }
    return s;
}

void ivanic_s_free(ivanic_s *s)
{
    if (s != NULL)
    {
        if (s->R != NULL)
        {
            my_matrix_free(s->R);
        }
        if (s->Rz_alpha != NULL)
        {
            my_matrix_free(s->Rz_alpha);
        }
        if (s->Ry_beta != NULL)
        {
            my_matrix_free(s->Ry_beta);
        }
        if (s->Rz_gamma != NULL)
        {
            my_matrix_free(s->Rz_gamma);
        }
        if (s->temp != NULL)
        {
            my_matrix_free(s->temp);
        }
    }
}

static double P(int i, int l, int mu, int m_prime, gsl_matrix *R, int R_lm1_offset)
{
    i += 2; // shift i by 2 to skip S0
    mu += l - 1 + R_lm1_offset;
    int m_prime_lm1 = m_prime + l - 1 + R_lm1_offset;
    int twolm2 = 2 * l - 2 + R_lm1_offset;
    if (abs(m_prime) < l)
    {
        return (double)my_matrix_get(R, i, 2) * (double)my_matrix_get(R, mu, m_prime_lm1);
    }
    if (m_prime == l)
    {
        return (double)my_matrix_get(R, i, 3) * (double)my_matrix_get(R, mu, twolm2) - (double)my_matrix_get(R, i, 1 + 0) * (double)my_matrix_get(R, mu, R_lm1_offset);
    }
    return (double)my_matrix_get(R, i, 3) * (double)my_matrix_get(R, mu, R_lm1_offset) + (double)my_matrix_get(R, i, 1 + 0) * (double)my_matrix_get(R, mu, twolm2);
}

static double U(int l, int m1, int m2, gsl_matrix *R, int R_lm1_offset)
{
    return P(0, l, m1, m2, R, R_lm1_offset);
}

static double V(int l, int m1, int m2, gsl_matrix *R, int R_lm1_offset)
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

static double W(int l, int m1, int m2, gsl_matrix *R, int R_lm1_offset)
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

static void
uvw(double buffer[3], int l, int m1, int m2)
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

    // Ry_beta
    fill_Ry(s->Ry_beta, cos_beta, sin_beta);

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
                    // compute_function_value(U, uvw_buffer[0], n, m1, m2, n_start, R, &M, P, mode);
                }
                // ------- V ------
                if (fabs(uvw_buffer[1]) > 1e-15)
                {
                    // compute_function_value(V, uvw_buffer[1], n, m1, m2, n_start, R, &M, P, mode);
                }
                // ------- W ------
                if (fabs(uvw_buffer[2]) > 1e-15)
                {
                    // compute_function_value(W, uvw_buffer[2], n, m1, m2, n_start, R, &M, P, mode);
                }
            }
        }
    }
}

typedef struct _mtx_spherical_harmonics_rotation_
{
    t_object x_obj;
    t_outlet *list_rot_out;
    t_atom *list_rot;
    ivanic_s *ivanic;
    double alpha;
    double beta;
    double gamma;
} t_mtx_spherical_harmonics_rotation;

static t_class *mtx_spherical_harmonics_rotation_class;

static t_mtx_spherical_harmonics_rotation *mtx_spherical_harmonics_rotation_new(t_symbol *s, int argc, t_atom *argv)
{
    post("[mtx_spherical_harmonics_rotation]: new called");
    int N;
    t_mtx_spherical_harmonics_rotation *x = (t_mtx_spherical_harmonics_rotation *)pd_new(mtx_spherical_harmonics_rotation_class);
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
        post("[mtx_spherical_harmonics_rotation]: n = %d", N);
    }
    x->list_rot = (t_atom *)calloc((size_t)pow((N + 1), 4) + 2, sizeof(t_atom));
    x->alpha = 0.0;
    x->beta = 0.0;
    x->gamma = 0.0;
    return (void *)x;
}

static void mtx_spherical_harmonics_rotation_free(t_mtx_spherical_harmonics_rotation *x)
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

static void mtx_spherical_harmonics_rotation_bang(t_mtx_spherical_harmonics_rotation *x)
{
    post("[mtx_spherical_harmonics_rotation]: bang called");
    if (x->list_rot != 0)
    {
        outlet_anything(x->list_rot_out, gensym("matrix"),
                        (size_t)pow((x->ivanic->N + 1), 4) + 2, x->list_rot);
    }
}

static void mtx_spherical_harmonics_rotation_list(t_mtx_spherical_harmonics_rotation *x, t_symbol *s,
                                                  int argc, t_atom *argv)
{
    size_t num_harmonics, n;
    post("mtx_spherical_harmonics_rotation_list called with %d args", argc);
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
    post("Angles: alpha %f beta %f gamma %f", x->alpha, x->beta, x->gamma);
    ivanic_ruedenberg_rotation_matrix(x->ivanic, x->alpha, x->beta, x->gamma);
    num_harmonics = (size_t)((x->ivanic->N + 1) * (x->ivanic->N + 1));
    SETFLOAT(x->list_rot, (t_float)num_harmonics);
    SETFLOAT(x->list_rot + 1, (t_float)num_harmonics);
    for (n = 0; n < num_harmonics * num_harmonics; ++n)
    {
        SETFLOAT(x->list_rot + n + 2, (t_float)x->ivanic->R->data[n]);
    }
    mtx_spherical_harmonics_rotation_bang(x);
}

/* this is called once at setup time, when this code is loaded into Pd.
 */
void mtx_spherical_harmonics_rotation_setup(void)
{
    mtx_spherical_harmonics_rotation_class = class_new(
        gensym("mtx_spherical_harmonics_rotation"),
        (t_newmethod)mtx_spherical_harmonics_rotation_new,
        (t_method)mtx_spherical_harmonics_rotation_free,
        sizeof(t_mtx_spherical_harmonics_rotation),
        CLASS_DEFAULT,
        A_GIMME, 0);
    class_addbang(mtx_spherical_harmonics_rotation_class, (t_method)mtx_spherical_harmonics_rotation_bang);
    class_addlist(mtx_spherical_harmonics_rotation_class, (t_method)mtx_spherical_harmonics_rotation_list);

#ifdef HAVE_LIBGSL
    my_matrix_calloc = iemmatrix_get_stub("gsl_matrix_calloc", mtx_spherical_harmonics_rotation_class);
    my_matrix_free = iemmatrix_get_stub("gsl_matrix_free", mtx_spherical_harmonics_rotation_class);
    my_matrix_get = iemmatrix_get_stub("gsl_matrix_get", mtx_spherical_harmonics_rotation_class);
    my_matrix_set = iemmatrix_get_stub("gsl_matrix_set", mtx_spherical_harmonics_rotation_class);
#endif
    post("[mtx_spherical_harmonics_rotation]: gsl stubs: calloc=%p free=%p get=%p set=%p",
         my_matrix_calloc, my_matrix_free, my_matrix_get, my_matrix_set);
    have_gsl = (my_matrix_calloc && my_matrix_free && my_matrix_get && my_matrix_set);
}
