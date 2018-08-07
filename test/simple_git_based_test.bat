cd ../bin
REM echo off
echo Testing primer3_core.exe

echo Test 1 : check '-output' and 'input_file'
primer3_core.exe  --default_version=1 --strict --output=../test/primer_check_output ../test/primer_check_input
echo Use git to see how the file primer_check_output changed
pause

echo Test 2 : check ' output' and ' input_file'
primer3_core.exe  --default_version=1 --st > ../test/primer_check_output < ../test/primer_check_input
echo Use git to see how the file primer_check_output changed
pause