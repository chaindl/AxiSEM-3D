//
//  Domain.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  computational domain

#ifndef Domain_hpp
#define Domain_hpp

#include <vector>
#include <memory>

// mesh
class Point;
class Element;

// boundary
class SolidFluidBoundary;
class AbsorbingBoundary;
class AxialBoundary;
class FluidSurfaceBoundary;

// mpi
class Messaging;

// source
class ElementSource;

// station output
class StationSolid;
class StationFluid;
template <class StationT>
class StationGroup;
typedef StationGroup<StationSolid> StationGroupInSolid;
typedef StationGroup<StationFluid> StationGroupInFluid;

// element output
class ElementOpSolid;
class ElementOpFluid;
template <class ElementT>
class ElementOpGroup;
typedef ElementOpGroup<ElementOpSolid> ElementOpGroupInSolid;
typedef ElementOpGroup<ElementOpFluid> ElementOpGroupInFluid;

// wavefield scanning
class WavefieldScanning;

class Domain {
public:
    // constructor
    Domain();
    
    ////////////////////// domain construction //////////////////////
    // add a point
    void addPoint(const std::shared_ptr<Point> &point);
    
    // add an element
    void addElement(const std::shared_ptr<Element> &element);
    
    // replace an element
    void replaceElement(const std::shared_ptr<Element> &element);
    
    // set mpi messaging
    void setMessaging(std::unique_ptr<Messaging> &msg);
    
    // verbose
    std::string verbose(const std::string &title) const;
    
    // MESHING STAGE //
    ///////////////////
    // SOLVING STAGE //
    
    // set wavefield injection
    // void setWavefieldInjection(std::unique_ptr<WavefieldInjection> &wj);
    
    // add a source
    void addElementSource(std::unique_ptr<const ElementSource> &esrc);
    
    // add station group in solid
    void addStationGroupInSolid(std::unique_ptr<StationGroupInSolid> &stgrp);
    
    // add station group in fluid
    void addStationGroupInFluid(std::unique_ptr<StationGroupInFluid> &stgrp);
    
    // add element group in solid
    void
    addElementOpGroupInSolid(std::unique_ptr<ElementOpGroupInSolid> &elgrp);
    
    // add element group in fluid
    void
    addElementOpGroupInFluid(std::unique_ptr<ElementOpGroupInFluid> &elgrp);
    
    // set wavefield scanning
    void setWavefieldScanning(std::unique_ptr<WavefieldScanning> &ws);
    
    
    ////////////////////// get domain //////////////////////
    // get solid points
    inline const std::vector<std::shared_ptr<Point>> &
    getPoints() const {
        return mPoints;
    }

    // get elements
    inline const std::vector<std::shared_ptr<const Element>> &
    getElements() const {
        return mElements;
    }
    
    // get solid-fluid boundary
    inline std::unique_ptr<SolidFluidBoundary> &getSolidFluidBoundary() {
        return mSolidFluidBoundary;
    }
    
    // get absorbing boundary
    inline std::unique_ptr<AbsorbingBoundary> &getAbsorbingBoundary() {
        return mAbsorbingBoundary;
    }
    
    // get axial boundary
    inline std::unique_ptr<AxialBoundary> &getAxialBoundary() {
        return mAxialBoundary;
    }
    
    // get fluid surface boundary
    inline std::unique_ptr<FluidSurfaceBoundary> &getFluidSurfaceBoundary() {
        return mFluidSurfaceBoundary;
    }
    
    ////////////////////// timeloop //////////////////////
    // apply sources
    void applySources(int tstep, double time) const;
    
    // compute stiffness
    void computeStiffness() const;
    
    // boundary conditions before assembling stiffness
    // * solid-fluild
    void applyBC_BeforeAssemblingStiff() const;
    
    // boundary conditions after assembling stiffness
    // * Clayton ABC
    // * axial
    // * fluid surface
    void applyBC_AfterAssemblingStiff() const;
    
    // boundary conditions after computing acceleration
    // * sponge ABC
    void applyBC_AfterComputingAccel() const;
    
    // stiff to accel on points
    void computeStiffToAccel() const;
    void combinePointWindows() const;
    void separatePointWindows() const;
    
    // mpi phase 1: commGatherSendRecv
    void mpiGatherSendRecv() const;
    
    // mpi phase 2: commWaitScatter
    void mpiWaitScatter() const;
    
    // initialize output
    void initializeOutput() const;
    
    // record output
    void recordOutput(int tstep, double time) const;
    
    // dump output
    void dumpOutput() const;
    
    // finalize output
    void finalizeOutput() const;
    
    // do scanning
    void doScanning(int tstep) const;
    
    // report scanning
    void reportScanning() const;
    
    // check stability
    void checkStability(int tstep, double t, double dt) const;
    
    
private:
    ////////////////////// spectral elements //////////////////////
    // points
    std::vector<std::shared_ptr<Point>> mPoints;
    
    // elements
    std::vector<std::shared_ptr<const Element>> mElements;
    
    
    ////////////////////// boundaries //////////////////////
    // solid-fluid boundary
    std::unique_ptr<SolidFluidBoundary> mSolidFluidBoundary;
    
    // absorbing boundary
    std::unique_ptr<AbsorbingBoundary> mAbsorbingBoundary;
    
    // axial boundary
    std::unique_ptr<AxialBoundary> mAxialBoundary;
    
    // fluid surface
    std::unique_ptr<FluidSurfaceBoundary> mFluidSurfaceBoundary;
    
    
    ////////////////////// mpi //////////////////////
    std::unique_ptr<Messaging> mMessaging;
    
    // MESHING STAGE //
    ///////////////////
    // SOLVING STAGE //
    
    ////////////////////// source //////////////////////
    // injection
    // std::unique_ptr<const WavefieldInjection> mWavefieldInjection;
    
    // sources
    std::vector<std::unique_ptr<const ElementSource>> mElementSources;
    
    
    ////////////////////// output //////////////////////
    // stations in solid
    std::vector<std::unique_ptr<StationGroupInSolid>> mStationGroupInSolids;
    
    // stations in fluid
    std::vector<std::unique_ptr<StationGroupInFluid>> mStationGroupInFluids;
    
    // elements in solid
    std::vector<std::unique_ptr<ElementOpGroupInSolid>> mElementOpGroupInSolids;
    
    // elements in fluid
    std::vector<std::unique_ptr<ElementOpGroupInFluid>> mElementOpGroupInFluids;
    
    
    ////////////////////// wavefield scanning //////////////////////
    std::unique_ptr<const WavefieldScanning> mWavefieldScanning;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // empty creator
    // constructor not available from outside becuase of forward declarations
    static std::shared_ptr<Domain> createEmpty();
};

#endif /* Domain_hpp */
