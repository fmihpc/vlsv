/** This file is part of VLSV file format. The file contains wrapper functions 
* for some of the functions defined in Linux header unistd.h. In Windows 
* systems define preprocessor macro "WINDOWS".
*
*  Copyright 2015 Arto Sandroos
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU Lesser General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#ifndef PORTABLE_FILE_SYSTEM_H
#define PORTABLE_FILE_SYSTEM_H

namespace fileio {

	int chdir(const char* path);
	char* getcwd(char* buf, size_t size);

} // namespace fileio

#endif