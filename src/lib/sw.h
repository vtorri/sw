#ifndef SW_H
#define SW_H


#ifdef SAPI
# undef SAPI
#endif

#ifdef _WIN32
# ifdef SW_BUILD
#  ifdef DLL_EXPORT
#   define SAPI __declspec(dllexport)
#  else
#   define SAPI
#  endif
# else
#  define SAPI __declspec(dllimport)
# endif
#else
# ifdef __GNUC__
#  if __GNUC__ >= 4
#   define SAPI __attribute__ ((visibility("default")))
#  else
#   define SAPI
#  endif
# else
#  define SAPI
# endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <sw_spline.h>
#include <sw_scale_fct.h>
#include <sw_weights.h>
#include <sw_sweldens.h>
#include <sw_wavelet.h>
#include <sw_mra.h>
#include <sw_utils.h>

#ifdef __cplusplus
}
#endif


#endif /* SW_H */
