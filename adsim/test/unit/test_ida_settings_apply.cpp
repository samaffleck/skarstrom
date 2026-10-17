// ApplyIDASettings must hand IDA only values it accepts, and report the first
// rejected call so the solver can fail instead of silently running on IDA's
// defaults (L-14).

#include <gtest/gtest.h>

#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>

#include "adsim/solver/ida_settings.h"

namespace adsim {

  namespace {

    // Trivial residual F = y' + y, only needed so IDAInit accepts the problem.
    int Residual(sunrealtype, N_Vector y, N_Vector yp, N_Vector r, void*) {
      N_VLinearSum(1.0, yp, 1.0, y, r);
      return 0;
    }

    // A real, initialised IDA instance (some setters, e.g. IDASetMaxOrd and
    // IDASetMaxNonlinIters, depend on state set up by IDAInit).
    class IdaInstance {
    public:
      IdaInstance() {
        SUNContext_Create(SUN_COMM_NULL, &ctx_);
        y_ = N_VNew_Serial(1, ctx_);
        yp_ = N_VNew_Serial(1, ctx_);
        N_VConst(1.0, y_);
        N_VConst(-1.0, yp_);
        mem_ = IDACreate(ctx_);
        init_flag_ = IDAInit(mem_, Residual, 0.0, y_, yp_);
      }
      ~IdaInstance() {
        IDAFree(&mem_);
        N_VDestroy(y_);
        N_VDestroy(yp_);
        SUNContext_Free(&ctx_);
      }
      void* mem() const { return mem_; }
      int init_flag() const { return init_flag_; }

    private:
      SUNContext ctx_ = nullptr;
      N_Vector y_ = nullptr;
      N_Vector yp_ = nullptr;
      void* mem_ = nullptr;
      int init_flag_ = IDA_MEM_NULL;
    };

  } // namespace

  TEST(IDASettingsApply, DefaultSettingsApplyCleanly) {
    IdaInstance ida;
    ASSERT_EQ(ida.init_flag(), IDA_SUCCESS);
    IDASettings settings;
    EXPECT_EQ(ApplyIDASettings(ida.mem(), settings), IDA_SUCCESS);
  }

  TEST(IDASettingsApply, ExplicitMaxStepAppliesCleanly) {
    IdaInstance ida;
    ASSERT_EQ(ida.init_flag(), IDA_SUCCESS);
    IDASettings settings;
    settings.max_step = 2.0;
    settings.max_order = 3;
    EXPECT_EQ(ApplyIDASettings(ida.mem(), settings), IDA_SUCCESS);
  }

  TEST(IDASettingsApply, RejectedSettingIsReported) {
    IdaInstance ida;
    ASSERT_EQ(ida.init_flag(), IDA_SUCCESS);
    IDASettings settings;
    settings.max_order = 6; // IDA allocates for at most BDF order 5
    EXPECT_EQ(ApplyIDASettings(ida.mem(), settings), IDA_ILL_INPUT);
  }

} // namespace adsim
