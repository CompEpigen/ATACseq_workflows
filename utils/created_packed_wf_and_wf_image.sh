DIR_OF_THIS_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${DIR_OF_THIS_SCRIPT}/../CWL/workflows/"
if [[ -d "./packed" ]]
then
    rm -r "./packed"
fi
if [[ -d "./graph_images" ]]
then
    rm -r "./graph_images"
fi
mkdir "./packed" "./graph_images"
CWL_FILES=( $( ls *.cwl ) )
for CWL_FILE in ${CWL_FILES[@]}
do
    cwltool --print-dot "$CWL_FILE"| dot -Tpdf > "./graph_images/${CWL_FILE}.pdf"
    cwltool --pack "$CWL_FILE" > "./packed/${CWL_FILE}"
done