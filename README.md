Medical rota design/checking tool.
By Rudolf Cardinal (rudolf@pobox.com), Aug 2015.

Copyright/licensing

    Copyright (C) 2015-2015 Rudolf Cardinal (rudolf@pobox.com).
    Department of Psychiatry, University of Cambridge.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

Purpose

    Fed up with financially correct but clinically poor rotas emanating from
    your Medical Staffing department?

    Able to write clinically safe rotas but unable to check them for validity
    and banding calculations?

    Perhaps this can help.

    The tool allows you to write rotas (at present, by adding a new function
    to the code, which is written in Python 3).

    You define:
        - shifts
        - doctors, each with a shift pattern (typically all doctors have a
          rotated copy of the same pattern)
        - a few rota-wide settings

    The software generates and displays the full rota for a given time period,
    and attempts to check

    - role coverage (do you have enough people on at night? Do you always have
      a junior doctor on call who's approved under section 12 of the Mental
      Health Act?)

    - coverage of the normal working day (which is the thing that suffers if
      you spend too much time on call)

    - banding, and New Deal/European Working Time Directive compliance

No guarantees

    I don't guarantee the software, but make it freely available in the hope
    that it saves you a few iterations begging your Medical Staffing department
    to check your various rota designs for banding/compliance, and thus allows
    you to develop something safe and financially reasonable in a shorter time
    and with fewer false starts.

Tips

    - The program produces HTML output; use a web browser to view the rota
      files generated.

    - Use the '-n' or '--daynums' option to add day numbers (as well as dates)
      to the rota. This can help when comparing existing/new rotas.

    - Use the '-c' or '--nocalc' option to generate a rota very quickly
      but without the analysis. Iterate through (change code, regenerate,
      check) quickly until you're satisfied; then remove this option to perform
      the full analyses.
