#ifndef AVT_VLSV_DATABASE_H
#define AVT_VLSV_DATABASE_H

#include <avtDatasetDatabase.h>
#include <avtGenericDatabase.h>
#include <avtMaterial.h>
#include <avtSpecies.h>
#include <avtVariableCache.h>
#include <avtTransformManager.h>
#include <avtFileFormatInterface.h>

class DATABASE_API avtVlsvDatabase : public avtGenericDatabase
{
    public:
        avtVlsvDatabase(avtFileFormatInterface *);
        virtual ~avtVlsvDatabase();
    protected:
        virtual void SetDatabaseMetaData(avtDatabaseMetaData *md, int timeState, bool forceReadAllCyclesTimes = true);
};

#endif