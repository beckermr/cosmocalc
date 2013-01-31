
#ifndef _COSMOCALC_ASSERT_
#define _COSMOCALC_ASSERT_
/* Some debugging macros
   undef DEBUG for no debugging
   DEBUG_LEVEL = 0 is for basic debugging
   DEBUG_LEVEL = 1 is for messages printed by a single task but not as critical
   DEBUG_LEVEL = 2 or above is used for messages that every task will print
*/

#ifdef NDEBUG
#undef DEBUG
#define DEBUG_LEVEL -1
#undef DEBUG_IO
#endif

#ifdef DEBUG
#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif
#endif

#ifdef NDEBUG
#define cosmocalc_assert(ex, format, ...) ((void)0)
#else
#include <cstdio>
#include <cstdlib>
#define cosmocalc_assert(ex, format, ...) do { \
    if(!(ex)) {\
      printf("cosmoCalc ERROR: %s:%u: '%s' - "format"\n", __FILE__, __LINE__, #ex, ## __VA_ARGS__); \
      exit(EXIT_FAILURE); \
    } \
  } while(false)
#endif

#endif /* _COSMOCALC_ASSERT_ */
