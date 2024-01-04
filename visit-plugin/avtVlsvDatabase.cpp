#include <avtVlsvDatabase.h>

#include <algorithm>
#include <float.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCSGGrid.h>
#include <vtkDataSet.h>
#include <vtkDataSetWriter.h>
#include <vtkDoubleArray.h>
#include <vtkEnumThreshold.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataRelevantPointsFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkTrivialProducer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridBoundaryFilter.h>
#include <vtkUnstructuredGridFacelistFilter.h>
#include <vtkVisItUtility.h>

#include <visitstream.h>

#include <avtCallback.h>
#include <avtDatabaseMetaData.h>
#include <avtDatasetCollection.h>
#include <avtDatasetVerifier.h>
#include <avtDebugDumpOptions.h>
#include <avtDomainBoundaries.h>
#include <avtDomainNesting.h>
#include <avtFileFormatInterface.h>
#include <avtGhostNodeGenerator.h>
#include <avtMemory.h>
#include <avtMixedVariable.h>
#include <avtParallel.h>
#include <avtSILGenerator.h>
#include <avtSILRestrictionTraverser.h>
#include <avtSourceFromDatabase.h>
#include <avtStreamingGhostGenerator.h>
#include <avtStructuredDomainBoundaries.h>
#include <avtLocalStructuredDomainBoundaries.h>
#include <avtStructuredDomainNesting.h>
#include <avtTransformManager.h>
#include <avtTypes.h>
#include <avtUnstructuredPointBoundaries.h>
#include <PickAttributes.h>
#include <PickVarInfo.h>
#include <QueryOverTimeAttributes.h>

#include <DebugStream.h>
#include <BadDomainException.h>
#include <ImproperUseException.h>
#include <InvalidDBTypeException.h>
#include <InvalidVariableException.h>
#include <NoInputException.h>
#include <TimingsManager.h>



using     std::string;
using     std::vector;

avtVlsvDatabase::avtVlsvDatabase(avtFileFormatInterface *inter)
: avtGenericDatabase(inter)
{
}

avtVlsvDatabase::~avtVlsvDatabase()
{
}

// ****************************************************************************
//  Method: avtGenericDatabase::SetDatabaseMetaData
//
//  Purpose:
//      Sets the database meta-data using the file format interface.
//
//  Arguments:
//      md        : The meta-data to set.
//      timeState : The time state that we're interested in.
//
//  Programmer: Hank Childs
//  Creation:   March 2, 2001
//
//  Modifications:
//    Brad Whitlock, Wed May 14 09:15:18 PDT 2003
//    I added the optional timeState argument.
//
//    Hank Childs, Tue Feb 15 07:21:10 PST 2005
//    Replace forbidden characters of expression language.
//
//    Mark C. Miller, Tue May 17 18:48:38 PDT 2005
//    Added bool arg, forceReadAllCyclesTimes
//
//    Jeremy Meredith, Fri Sep  2 15:03:06 PDT 2005
//    Removed most of the special character replacements, as the expression
//    language scanner is now accepting most of them inside <>'s.  Replaced
//    the unacceptable [] <> () with the acceptable {}'s since it seems
//    more natural than the former _ character used to replace them.
//
//    Hank Childs, Sun Mar 16 07:10:49 PDT 2008
//    Change substitutions for (, }, [, ], <, and > to use an underscore,
//    instead of { & }, since that was screwing up expressions.
//
//    Dave Pugmire, Thu Mar  4 10:27:17 EST 2010
//    Ensure that naming collisions don't occur. Previously ()[]<> all mapped
//    to '_'.
//
//    Markku Alho: Force reading time & cycle from vlsv files
// ****************************************************************************

void avtVlsvDatabase::SetDatabaseMetaData(avtDatabaseMetaData *md, int timeState, bool forceReadAllCyclesTimes)
{
    debug2 << "VLSV\t SetDataBaseMetaData called\n";
    int t0 = visitTimer->StartTimer();
    Interface->SetDatabaseMetaData(md, timeState, true);
    md->ReplaceForbiddenCharacters();
    visitTimer->StopTimer(t0, "Getting database meta data");
}