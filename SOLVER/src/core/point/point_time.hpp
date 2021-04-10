//
//  point_time.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/11/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  time-scheme-dependent point methods

#ifndef point_time_hpp
#define point_time_hpp

class Point;
class PointWindow;
class TimeScheme;

namespace point_time {
    // create fields
    void createFields(PointWindow &pw, const TimeScheme &timeScheme);
    
    // measure cost of a point
    double measure(Point &point, int count, const TimeScheme &timeScheme);
}

#endif /* point_time_hpp */
