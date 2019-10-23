#ifndef CNTR_RESPONSE_CONVOLUTION_IMPL_H
#define CNTR_RESPONSE_CONVOLUTION_IMPL_H

#include "cntr_response_convolution_decl.hpp"

namespace cntr {

/*////////////////////////////////////////////////////////////////

        computing the response to a retarded response function

        c(t) =  \int_C  dt' W(t,t') f(t')

        NOTES:
        - assumes that W is hermitian
        - uses points 0...kt if t<kt (thus a version using W only as timestep
          works only for t>=kt and t==-1

////////////////////////////////////////////////////////////////*/

#define CPLX std::complex<T>

// res = int_0^n dt f(t),   f(i*dt)=f[sf*i+idxf], i=0..,n
/// @private
template <typename T>
void response_integrate(int n, T dt, CPLX &res, CPLX *f, int sf, int idxf,
                        integration::Integrator<T> &I) {
    int kt = I.k(), i, k2 = 2 * (kt + 1);
    res = 0;
    if (n >= k2) {
        for (i = 0; i <= kt; i++)
            res += f[i * sf + idxf] * I.gregory_weights(n, i);
        for (i = kt + 1; i < n - kt; i++)
            res += f[i * sf + idxf];
        for (i = n - kt; i <= n; i++)
            res += f[i * sf + idxf] * I.gregory_weights(n, i);
    } else if (n > 0) {
        int n1 = (n <= kt ? kt : n);
        for (i = 0; i <= n1; i++)
            res += f[i * sf + idxf] * I.gregory_weights(n, i);
    }
    res *= dt;
}
// res = int_0^n dt f(t)g(t),   f(i*dt)=f[sf*i+idxf], g(i*dt)=g[sg*i+idxg], i=0..,n
/// @private
template <typename T>
void response_integrate(int n, T dt, CPLX &res, CPLX *f, int sf, int idxf, CPLX *g, int sg,
                        int idxg, integration::Integrator<T> &I) {
    int kt = I.k(), i, k2 = 2 * (kt + 1);
    res = 0;
    if (n >= k2) {
        for (i = 0; i <= kt; i++)
            res += f[i * sf + idxf] * g[i * sg + idxg] * I.gregory_weights(n, i);
        for (i = kt + 1; i < n - kt; i++)
            res += f[i * sf + idxf] * g[i * sg + idxg];
        for (i = n - kt; i <= n; i++)
            res += f[i * sf + idxf] * g[i * sg + idxg] * I.gregory_weights(n, i);
    } else if (n > 0) {
        int n1 = (n <= kt ? kt : n);
        for (i = 0; i <= n1; i++)
            res += f[i * sf + idxf] * g[i * sg + idxg] * I.gregory_weights(n, i);
    }
    res *= dt;
}

// USING THE ASSUMPTION THAT W IS HERMITIAN
template <typename T>
void response_convolution(int tstp, CPLX &cc, GREEN_TSTP &W, int a1, int a2, function<T> &f,
                          int b1, int b2, int SolveOrder, T beta, T h) {
    int ntau = W.ntau();
    T dtau = beta / ntau;
    int sizew = W.size1(), idxw = a1 * sizew + a2, sw = sizew * sizew;
    int sizef = f.size1(), idxf = b1 * sizef + b2, sf = sizef * sizef;
    CPLX res;

    assert(W.tstp() >= tstp);
    assert(tstp == -1 || tstp >= SolveOrder);
    assert(SolveOrder <= ntau);
    assert(b1 >= 0 && b1 <= sizef-1);
    assert(b2 >= 0 && b2 <= sizef-1);
    assert(a1 >= 0 && a1 <= sizew-1);
    assert(a2 >= 0 && a2 <= sizew-1);

    if (tstp == -1) {
        response_integrate<T>(ntau, dtau, res, W.matptr(0), sw, idxw, integration::I<T>(SolveOrder));
        cc = res * f.ptr(-1)[idxf];
    } else {
        response_integrate<T>(ntau, dtau, res, W.tvptr(0), sw, idxw, integration::I<T>(SolveOrder));
        cc = res * (CPLX(0, -1.0) * f.ptr(-1)[idxf]);
        response_integrate<T>(tstp, h, res, W.retptr(0), sw, idxw, f.ptr(0), sf, idxf,
                              integration::I<T>(SolveOrder));
        cc += res;
    }
}
template <typename T>
void response_convolution(int tstp, CPLX &cc, GREEN &W, int a1, int a2, function<T> &f,
                          int b1, int b2, int SolveOrder, T beta, T h) {
    int ntau = W.ntau();
    T dtau = beta / ntau;
    int sizew = W.size1(), idxw = a1 * sizew + a2, sw = sizew * sizew;
    int sizef = f.size1(), idxf = b1 * sizef + b2, sf = sizef * sizef;
    int n1 = (tstp == -1 || tstp >= SolveOrder ? SolveOrder : SolveOrder);
    CPLX res;

    assert(W.nt() >= tstp);
    assert(tstp == -1 || tstp >= SolveOrder);
    assert(SolveOrder <= ntau);
    assert(b1 >= 0 && b1 <= sizef-1);
    assert(b2 >= 0 && b2 <= sizef-1);
    assert(a1 >= 0 && a1 <= sizew-1);
    assert(a2 >= 0 && a2 <= sizew-1);

    if (tstp == -1) {
        response_integrate<T>(ntau, dtau, res, W.matptr(0), sw, idxw, integration::I<T>(SolveOrder));
        cc = res * f.ptr(-1)[idxf];
    } else if (tstp >= SolveOrder) {
        response_integrate<T>(ntau, dtau, res, W.tvptr(tstp, 0), sw, idxw,
                              integration::I<T>(SolveOrder));
        cc = res * (CPLX(0, -1.0) * f.ptr(-1)[idxf]);
        response_integrate<T>(tstp, h, res, W.retptr(tstp, 0), sw, idxw, f.ptr(0), sf, idxf,
                              integration::I<T>(SolveOrder));
        cc += res;
    } else {
        response_integrate<T>(ntau, dtau, res, W.tvptr(tstp, 0), sw, idxw,
                              integration::I<T>(SolveOrder));
        cc = res * (CPLX(0, -1.0) * f.ptr(-1)[idxf]);
        // need to extrapolate W(t,t') for t'>t
        // ** USING THE ASSUMPTION THAT W IS HERMITIAN**
        // i.e., Wret(t,t') = -Wret(t',t)^*
        CPLX *wtmp = new CPLX[SolveOrder + 1];
        for (int i = 0; i <= tstp; i++)
            wtmp[i] = W.retptr(tstp, i)[a1 * sizew + a2];
        for (int i = tstp + 1; i <= SolveOrder; i++)
            wtmp[i] = -conj(W.retptr(i, tstp)[a2 * sizew + a1]);
        response_integrate<T>(tstp, h, res, wtmp, 1, 0, f.ptr(0), sf, idxf,
                              integration::I<T>(SolveOrder));
        cc += res;
        delete[] wtmp;
    }
}

#undef GREEN
#undef GREEN_TSTP
#undef CPLX

} // namespace cntr

#endif  // CNTR_RESPONSE_CONVOLUTION_IMPL_H
