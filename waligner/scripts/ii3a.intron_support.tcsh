#!bin/tcsh -f

set run=$1
set lane=$2
set pair=$3

set paired=""
if ($pair > 0) set  paired="-pair $pair"

echo "bin/bestali  -intronSupport -mrnaRemap tmp/METADATA/mrnaRemap.gz -i tmp/COUNT/$lane.hits.gz -o tmp/INTRONLANES/$lane  -run $run -gzo -unique $paired -minIntronOverlap $minIntronOverlap"
      bin/bestali  -intronSupport -mrnaRemap tmp/METADATA/mrnaRemap.gz -i tmp/COUNT/$lane.hits.gz -o tmp/INTRONLANES/$lane  -run $run -gzo -unique $paired -minIntronOverlap $minIntronOverlap
echo "bin/bestali  -intronSupport -mrnaRemap tmp/METADATA/mrnaRemap.gz -i tmp/COUNT/$lane.hits.gz -o tmp/INTRONLANES/$lane  -run $run -gzo $paired -minIntronOverlap $minIntronOverlap"
      bin/bestali  -intronSupport -mrnaRemap tmp/METADATA/mrnaRemap.gz -i tmp/COUNT/$lane.hits.gz -o tmp/INTRONLANES/$lane  -run $run -gzo $paired -minIntronOverlap $minIntronOverlap
echo 
exit 0
