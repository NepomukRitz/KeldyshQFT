#ifndef FPP_MFRG_RUNTIMEERROR_H
#define FPP_MFRG_RUNTIMEERROR_H

#include <stdexcept>
#include <string>

class RuntimeError : public std::runtime_error {
public:
    explicit RuntimeError(std::string const &reason);
    virtual ~RuntimeError() throw ();
private:
}; // RuntimeError

inline RuntimeError::RuntimeError(std::string const &reason)
        : std::runtime_error(reason) { }

inline RuntimeError::~RuntimeError() throw () { }



#endif //FPP_MFRG_RUNTIMEERROR_H
