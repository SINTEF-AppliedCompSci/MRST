#include <mex.h>

#include <sstream>
#include <system_error>

namespace {
    bool args_ok(const int nlhs, const int nrhs, const mxArray* prhs[])
    {
        return (nlhs <= 1) && (nrhs == 1)
            && mxIsScalar(prhs[0])
            && (mxIsDouble(prhs[0]) || mxIsInt32(prhs[0]));
    }

    int getExitCode(const mxArray* ecode)
    {
        return static_cast<int>(mxGetScalar(ecode));
    }

    std::string message(const int ecode)
    {
        return std::system_category().message(ecode);
    }

    mxArray* copyOut(std::string&& msg)
    {
        return mxCreateString(msg.c_str());
    }

    mxArray* translate(const mxArray* ecode)
    {
        return copyOut(message(getExitCode(ecode)));
    }
}

// Attempt to Translate a System Exit Code to Readable Message
//
// Syntax:
//   msg = mrstTranslateExitCode(ecode)

void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    if (args_ok(nlhs, nrhs, prhs)) {
        plhs[0] = translate(prhs[0]);
    }
    else {
        std::ostringstream os;
        os << "Syntax is\n\t"
           << "msg = " << mexFunctionName() << "(ecode);";

        mexErrMsgIdAndTxt("Input:Invalid", os.str().c_str());
    }
}
