#ifndef __statsAccumulator__
#define __statsAccumulator__

#include <string>

/*! \brief Base class for stats accumulators */
class statsAccumulator {
protected:
    /*! \brief Utility function for flux density setup */
    void setupFluxes(float, float, bool, unsigned int, float* s);

    /*! \brief Utility function for stats accumulation */
    void accumulateStats(unsigned int, unsigned int, 
                         float* const, float* const, 
                         float* const, float* const);

public:
    static const unsigned int nprob = 3;
    const float plevels[nprob] = {0.683, 0.954, 0.997};

    statsAccumulator() {}
    virtual ~statsAccumulator() {}

    /*! \breif Main processing function; constructs stats */
    virtual void build(const std::string&, unsigned int, bool) = 0;

    /*! \brief Serialize to file as HDF5 */
    virtual void writeAsHDF5(const std::string&) const = 0;
};

#endif