#include "cosmocalc.h"
#include "weaklens.h"

cosmocalcData cosmoData;
weaklensData wlData;

#include <gsl/gsl_errno.h>
void turn_off_gsl_errs(void) {
  gsl_set_error_handler_off();
}
