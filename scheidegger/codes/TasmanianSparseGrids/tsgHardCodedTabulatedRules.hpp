/*
 * Code Author: Miroslav Stoyanov
 *
 * Copyright (C) 2015  Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive Approximation
 *              a.k.a. TASMANIAN
 *
 * TASMANIAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TASMANIAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TASMANIAN.  If not, see <http://www.gnu.org/licenses/>
 *
 */

#ifndef __TASMANIAN_SPARSE_HARDCODED_RULES_HPP
#define __TASMANIAN_SPARSE_HARDCODED_RULES_HPP

#include "tsgCoreOneDimensional.hpp"

namespace TasGrid{

class TableGaussPatterson{
public:
        TableGaussPatterson();
        ~TableGaussPatterson();

        int getMaxLeve() const;
        double* getNodes( int level ) const;
        double getWeight( int level, int point ) const;

protected:
        void loadNodes();
        void loadWeights();

private:
        int max_levels;
        double *nodes; // contains the x-coordinate of each sample point
        double *weights; // contains the weight associated with each level

        int *weights_offsets;

        OneDimensionalMeta meta;
};

}

#endif
