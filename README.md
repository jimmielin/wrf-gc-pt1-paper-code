# WRF-GC | A one-way coupled meteo/chem model

This repository holds the `chem` directory for the WRF-GC model, which contains all the WRF-GC coupling code. **It supports GEOS-Chem 12.2.1**. The copy of GEOS-Chem is not included.

**This version of WRF-GC is for permanent archival for submission to Geosci. Model Dev.**
Please do not use this version for research. It may be outdated. Please refer to the official website at http://wrf.geos-chem.org to obtain the latest release.

(c) 2017-2019 Haipeng Lin <linhaipeng@pku.edu.cn> (<hplin@seas.harvard.edu>), Xu Feng <fengx7@pku.edu.cn>, Tzung-May Fu <fuzm@sustc.edu.cn>

GEOS-Chem, GEOS-Chem High Performance, HEMCO, ESMF/MAPL Frameworks are (c) their original authors.

## Primary Components for Coupling

* WRF-to-Chemistry Abstraction Layer (Codename "Pumpkin") - `wrf-gchp-pumpkin` Project
* WRF-GC Chemistry Driver - `chemics_init` & `chem_driver`
* Stateful Conversion Module - `WRFGC_Convert_State_Mod` (formerly `GIGC_Convert_State_Mod`)

## Support Components
The below support components are based off GEOS-Chem High Performance ("GCHP") technology.

* GEOS-Chem Column Code - `GIGC_Chunk_Mod`, based off GCHP's original column code
* GEOS-Chem Stateful Variables Container - `GC_Stateful_Mod`, designed by WRF-GC project to support multiple domains within the same CPU in GEOS-Chem

## License
```
(c) 2017-2019 Haipeng Lin <hplin@seas.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>, Tzung-May Fu <fuzm@sustc.edu.cn>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to 
use, copy, modify the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

- The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

- The Software, modified in part or in full may not be redistributed without
express permission from the copyright holder(s).

Except as contained in this notice or in attribution, the name of the WRF-GC model
shall not be used as an endorsement for distributing modified copies of the
Software without prior written permission from the copyright holder(s).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```