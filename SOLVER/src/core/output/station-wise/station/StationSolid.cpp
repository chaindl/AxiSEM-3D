//
//  StationSolid.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/3/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  station in solid

#include "StationSolid.hpp"
#include "geodesy.hpp"

/////////////////////////// setup ///////////////////////////
// set in group
void StationSolid::
setInGroup(int dumpIntv, const channel::solid::ChannelOptions &chops) {
    // member buffers
    if (chops.mNeedBufferU) {
        mBufferU.resize(dumpIntv, 3);
        mMajorityDisplInRTZ = mElement->getMajorityDisplInRTZ(mWindowPhis);
    }
    if (chops.mNeedBufferG) {
        mBufferG.resize(dumpIntv, 9);
        mMajorityNablaInRTZ = mElement->getMajorityNablaInRTZ(mWindowPhis);
    }
    if (chops.mNeedBufferE) {
        mBufferE.resize(dumpIntv, 6);
        mMajorityStrainInRTZ = mElement->getMajorityStrainInRTZ(mWindowPhis);
    }
    if (chops.mNeedBufferR) {
        mBufferR.resize(dumpIntv, 3);
        mMajorityCurlInRTZ = mElement->getMajorityCurlInRTZ(mWindowPhis);
    }
    if (chops.mNeedBufferS) {
        mBufferS.resize(dumpIntv, 6);
        mMajorityStressInRTZ = mElement->getMajorityStressInRTZ(mWindowPhis);
    }
    
    // element
    for (int m = 0; m < mWindowPhis.size(); m++) {
        mElement->prepareWavefieldOutput(chops, std::get<0>(mWindowPhis[m]), false);
    }
    
    // workspace
    int maxNu_1 = 0;
    for (int m = 0; m < mWindowPhis.size(); m++) {
        maxNu_1 = std::max(maxNu_1, mElement->getWindowNu_1(std::get<0>(mWindowPhis[m])));
    }
    expandWorkspaceRecord(maxNu_1, chops);
    expandWorkspaceProcess(dumpIntv, chops.mNeedBufferE || chops.mNeedBufferS);
}


/////////////////////////// record ///////////////////////////
// record
void StationSolid::
record(int bufferLine, const channel::solid::ChannelOptions &chops) {
    for (int m = 0; m < mWindowPhis.size(); m++) {
        int nu_1;
        
        // displ
        if (chops.mNeedBufferU) {
            mElement->getDisplField(sUXN3, mMajorityDisplInRTZ, std::get<0>(mWindowPhis[m]));
            interpolate<3>(sUXN3, sUX3, sU3, mElement->getWindowNu_1_noBuffer(std::get<0>(mWindowPhis[m])), m);
            mBufferU.row(bufferLine) += sU3;
        }
        // nabla
        if (chops.mNeedBufferG) {
            mElement->getNablaField(sGXN9, mMajorityNablaInRTZ, std::get<0>(mWindowPhis[m]), nu_1);
            interpolate<9>(sGXN9, sGX9, sG9, nu_1, m);
            mBufferG.row(bufferLine) += sG9;
        }
        // strain
        if (chops.mNeedBufferE) {
            mElement->getStrainField(sEXN6, mMajorityStrainInRTZ, std::get<0>(mWindowPhis[m]), nu_1);
            interpolate<6>(sEXN6, sEX6, sE6, nu_1, m);
            mBufferE.row(bufferLine) += sE6;
        }
        // curl
        if (chops.mNeedBufferR) {
            mElement->getCurlField(sRXN3, mMajorityCurlInRTZ, std::get<0>(mWindowPhis[m]), nu_1);
            interpolate<3>(sRXN3, sRX3, sR3, nu_1, m);
            mBufferR.row(bufferLine) += sR3;
        }
        // stress
        if (chops.mNeedBufferS) {
            mElement->getStressField(sSXN6, mMajorityStressInRTZ, std::get<0>(mWindowPhis[m]), nu_1);
            interpolate<6>(sSXN6, sSX6, sS6, nu_1, m);
            mBufferS.row(bufferLine) += sS6;
        }
    }
}


/////////////////////////// process ///////////////////////////
// process and report to group
void StationSolid::
processReport(int bufferLine,
              const channel::solid::ChannelOptions &chops,
              int stationIndex, eigen::RTensor3 &bufferFields) {
    // rotate
    rotate(bufferLine, chops);
    
    // loop over channels
    for (int ich = 0; ich < chops.mStdChannels.size(); ich++) {
        // find field and index of channel
        const int &cha = chops.mStdChannels[ich];
        const auto &tup = channel::solid::gChannelMap.at(cha);
        channel::solid::FieldType ftype = std::get<2>(tup);
        int fieldIndex = std::get<3>(tup);
        // compute and feed
        if (ftype == channel::solid::FieldType::Displ) {
            computeFeedChannel<3>(mBufferU, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::solid::FieldType::Nabla) {
            computeFeedChannel<9>(mBufferG, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::solid::FieldType::Strain) {
            computeFeedChannel<6>(mBufferE, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::solid::FieldType::Curl) {
            computeFeedChannel<3>(mBufferR, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::solid::FieldType::Stress) {
            computeFeedChannel<6>(mBufferS, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else {
            throw std::runtime_error("StationSolid::processReport || "
                                     "Unknown field type.");
        }
    }
}

// process 1: rotate
void StationSolid::
rotate(int bufferLine, const channel::solid::ChannelOptions &chops) {
    bool cartesian = geodesy::isCartesian();
    if (chops.mNeedBufferU) {
        rotateField<3>(mBufferU, bufferLine, mMajorityDisplInRTZ,
                       chops.mWCS, cartesian);
    }
    if (chops.mNeedBufferG) {
        rotateField<9>(mBufferG, bufferLine, mMajorityNablaInRTZ,
                       chops.mWCS, cartesian);
    }
    if (chops.mNeedBufferE) {
        // halve off-diagonal components before rotation
        mBufferE.rightCols(3) *= (numerical::Real).5;
        rotateField<6>(mBufferE, bufferLine, mMajorityStrainInRTZ,
                       chops.mWCS, cartesian);
        mBufferE.rightCols(3) *= (numerical::Real)2.;
    }
    if (chops.mNeedBufferR) {
        rotateField<3>(mBufferR, bufferLine, mMajorityCurlInRTZ,
                       chops.mWCS, cartesian);
    }
    if (chops.mNeedBufferS) {
        rotateField<6>(mBufferS, bufferLine, mMajorityStressInRTZ,
                       chops.mWCS, cartesian);
    }
}
