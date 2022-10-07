#include <iostream>
#include <algorithm>
#include "ReadManager.h"

using namespace std;
namespace bbpReader {

    Header::Header(char* readBuffer) {
#define fetch(TYPE, POSITION) *((TYPE*)(readBuffer + POSITION))
#define fetchS(POSITION) string((char*)(readBuffer + POSITION))

        double readIdentifier = fetch(double, IDENTIFIER_POSITION);
        headerSize = fetch(int32_t, HEADER_SIZE_POSITION);
        amountOfCells = fetch(int32_t, TOTAL_NUMBER_OF_CELLS_POSITION);
        /** \bug 32 bit is consistent with the binary format specification of the
            header, however
         */
        uint32_t compartmentsPerFrame = fetch(int32_t, TOTAL_NUMBER_OF_COMPARTMENTS_POSITION);
        libraryVersion = fetchS(LIBRARY_VERSION_POSITION);
        simulatorVersion = fetchS(SIMULATOR_VERSION_POSITION);
        numberOfSteps = fetch(int32_t, NUMBER_OF_STEPS_POSITION);
        startTime = fetch(double, TIME_START_POSITION);
        endTime = fetch(double, TIME_END_POSITION);
        dt = fetch(double, DT_TIME_POSITION);
        dataUnit = fetchS(D_UNIT_POSITION);
        timeUnit = fetchS(T_UNIT_POSITION);
        mappingSize = fetch(int32_t, MAPPING_SIZE_POSITION);
        mappingName = fetchS(MAPPING_NAME_POSITION);
        extraMappingSize = fetch(int32_t, EXTRA_MAPPING_SIZE_POSITION);
        extraMappingName = fetchS(EXTRA_MAPPING_NAME_POSITION);
        targetName = fetchS(TARGET_NAME_POSITION);

#undef fetchS
#undef fetch

        isNative = readIdentifier == ARCHITECTURE_IDENTIFIER;

        if (!isNative) {
            // check if data can be converted:
            DOUBLESWAP(readIdentifier);
            if (readIdentifier != ARCHITECTURE_IDENTIFIER) {
                cerr << "File is corrupt or originated from an unknown architecture." << endl;
                exit(1);
            }

            // if data was stored on a machine with different architecture,
            // then the data needs to be converted to be compatible with the reading machine
            FIX_INT(headerSize);
            FIX_INT(amountOfCells);
            FIX_INT(compartmentsPerFrame);
            FIX_INT(numberOfSteps);
            DOUBLESWAP(startTime);
            DOUBLESWAP(endTime);
            DOUBLESWAP(dt);
            FIX_INT(mappingSize);
            FIX_INT(extraMappingSize);
        }

        // printf( "read %d compartment count\n", compartmentsPerFrame );
        totalAmountOfCompartments = compartmentsPerFrame;
    }

    bool Header::operator==(const Header& h1) const {
#define _CHECK(field)      \
    if (h1.field != field) \
    return false
        // compare all fields other than isNative
        _CHECK(headerSize);
        _CHECK(libraryVersion);
        _CHECK(simulatorVersion);
        _CHECK(amountOfCells);
        _CHECK(totalAmountOfCompartments);
        _CHECK(numberOfSteps);
        _CHECK(startTime);
        _CHECK(endTime);
        _CHECK(dt);
        _CHECK(dataUnit);
        _CHECK(timeUnit);
        _CHECK(mappingSize);
        _CHECK(mappingName);
        _CHECK(extraMappingName);
        _CHECK(extraMappingSize);
        _CHECK(targetName);

        // else if nothing returned false:
        return true;
#undef _CHECK
    }

    CellInfo::CellInfo() {
        cellNum = 0;
        master = &cellNum;
        amountOfCompartments = 0;
        dataLocation = 0;
        extraMappingLocation = 0;
        mappingLocation = 0;
    }

    CellInfo::CellInfo(char* readBuffer, bool isNative, key cmp) {
#define fetch(TYPE, POSITION) *((TYPE*)(readBuffer + POSITION))
        cellNum = fetch(int32_t, NUMBER_OF_CELL_POSITION);
        amountOfCompartments = fetch(int32_t, NUMBER_OF_COMPARTMENTS_POSITION);

        dataLocation = fetch(uint64_t, DATA_INFO_POSITION);
        extraMappingLocation = fetch(uint64_t, EXTRA_MAPPING_INFO_POSITION);
        mappingLocation = fetch(uint64_t, MAPPING_INFO_POSITION);

#undef fetch

        if (!isNative) {
            // cout << "Fixing cell info data..." << endl;
            FIX_INT(cellNum);
            FIX_INT(amountOfCompartments);

            // cout << "size of offset type is " << sizeof(dataLocation) << endl;
            DOUBLESWAP(dataLocation);
            DOUBLESWAP(extraMappingLocation);
            DOUBLESWAP(mappingLocation);
        }

        // fprintf( stderr, "build cell info: gid %d dataPos %d extraPos %d mappingPos %d\n",
        // cellNum, (int) dataLocation, (int) extraMappingLocation, (int) mappingLocation );
        // cerr<<"build cell info: gid "<<cellNum<<" dataPos "<<dataLocation<<" extraPos
        // "<<extraMappingLocation<<" mappingPos "<<mappingLocation<<"\n";

        if (cmp == GID)
            master = &cellNum;
        else
            master = &dataLocation;
    }

    bool CellInfo::operator==(const CellInfo& c1) const {
#define _CHECK(field)      \
    if (c1.field != field) \
    return false  //; cout << # field << " is ok." << endl
        _CHECK(cellNum);
        _CHECK(amountOfCompartments);
        //_CHECK(frameSize);
        //_CHECK(dataSize);
        //_CHECK(extraMappingSize);
        _CHECK(dataLocation);
        _CHECK(extraMappingLocation);
        _CHECK(mappingLocation);

        return true;
#undef _CHECK
    }

    bool CellInfo::operator!=(const CellInfo& c1) const {
        return !(*this == c1);
    }

#define COMPARE_(op) (master == &cellNum) ? cellNum op c.cellNum : dataLocation op c.dataLocation

#define DEFINE_CMP_OP(op)                                 \
    bool CellInfo::operator op(const CellInfo& c) const { \
        return COMPARE_(op);                              \
    }

    DEFINE_CMP_OP(<)
    DEFINE_CMP_OP(<=)
    DEFINE_CMP_OP(>)
    DEFINE_CMP_OP(>=)

#undef DEFINE_CMP_OP
#undef COMPARE_

    bool CellInfo::compareByOffset(const CellInfo& a, const CellInfo& b) {
        return a.dataLocation < b.dataLocation;
    }

    FrameInfo::FrameInfo(char* buffer, char* endOfBuffer, bool isNative, CellInfo::key sortMode) {
        // create a cellinfo and then insert into itself (rem, this class inherits from set).  The
        // resulting set iterator is then put into
        // a map member for associative lookup gid:location
        for (char* bufferLoc = buffer; bufferLoc < endOfBuffer;
             bufferLoc += SIZE_CELL_INFO_LENGTH) {
            iterator ins = insert(CellInfo(bufferLoc, isNative, sortMode)).first;
            gidLookup.insert(make_pair(ins->getCellNum(), ins));
        }
    }

    FrameInfo::FrameInfo(const FrameInfo& rhs) : baseTypes::FrameInfo() {
        *this = rhs;
    }

    FrameInfo& FrameInfo::operator=(const FrameInfo& rhs) {
        insert(rhs.begin(), rhs.end());
        for (iterator cellItem = begin(); cellItem != end(); ++cellItem) {
            gidLookup.insert(make_pair(cellItem->getCellNum(), cellItem));
        }

        return *this;
    }

    FrameInfo& FrameInfo::filterGIDs(const vector<CellID>& gids) {
        Lookup newGidLookup;

        for (vector<CellID>::const_iterator gid = gids.begin(); gid != gids.end(); gid++) {
            Lookup::iterator it = gidLookup.find(*gid);
            if (it != gidLookup.end()) {
                newGidLookup.insert(*it);
                gidLookup.erase(it);
            } else
                throw string(
                    "GID doesn't exist in file.");  // TODO: replace by custom std::exception
        }

        // gidLookup contains the gids that need to be removed and the newGidLookup contains the
        // gids that remain.

        for (Lookup::iterator i = gidLookup.begin(); i != gidLookup.end(); i++)
            erase(i->second);  // delete the unwanted gids.

        gidLookup.swap(newGidLookup);

        return *this;
    }

    CompartmentMapping::CompartmentMapping() : Base() {
    }

    CompartmentMapping::CompartmentMapping(MappingItem* start, MappingItem* _end)
        : Base(start, _end) {
    }

    CompartmentMapping::CompartmentMapping(MappingItem* buffer, SizeType bufSize)
        : Base(buffer, buffer + bufSize) {
    }

    SectionData::SectionData() : Base() {
    }

    SectionData::SectionData(DataItem firstData) {
        push_back(firstData);
    }

    SectionData::SectionData(DataItem* start, DataItem* _end) : Base(start, _end) {
    }

    bool SectionData::operator==(const SectionData& s1) {
        if (s1.size() != size())
            return false;

        const_iterator dat1, dat2;
        for (dat1 = s1.begin(), dat2 = begin(); dat2 != end(); dat1++, dat2++)
            if (*dat1 != *dat2)
                return false;

        return true;
    }

    bool SectionData::operator!=(const SectionData& s1) {
        return !(*this == s1);
    }

    bool SectionData::compareByOriginalLocation(const Identifier& c1, const Identifier& c2) {
        return c1.getLocation() < c2.getLocation();
    }

    CellData::CellData() {
        sectSort = SID;
    }

    CellData::CellData(SortingKey sorting) {
        sectSort = sorting;
    }

    CellData::CellData(const DataItem* dataBuffer, const CompartmentMapping& mapping) {
        sectSort = SID;
        fillData(dataBuffer, mapping);
    }

    void CellData::fillData(const DataItem* dataBuffer, const CompartmentMapping& mapping) {
        CompartmentMapping::const_iterator it;
        SectionIndex i;
        for (it = mapping.begin(), i = 0; it != mapping.end(); it++, dataBuffer++, i++) {
            // iteration through all compartments inside cell
            DataItem data = *dataBuffer;

            iterator segMatch = find(*it);  // find corresponding segment

            if (segMatch == end()) {
                // SID was not found.
                // this means it wasnt created yet.
                // therefore create a new one and insert first data item into it:
                iterator ins =
                    insert(make_pair(SectionData::Identifier(*it, i, sectSort), SectionData(data)))
                        .first;
                // note that it is guaranteed that the above insertion takes place, because this
                // block is only entered
                // if the sid does not exist inside this CellData object yet.
                // if sorting is done by location, this will be unique anyways.
                sidLookup.insert(make_pair(*it, ins));
            } else
                // SID is already present in CellData. The new data item just needs to be inserted
                segMatch->second.push_back(data);
        }
    }

    // typedef CellDataImpl<SectionData::LOCATION> UnsortedCellData;

    bool CellData::operator==(const CellData& d1) {
        if (d1.size() != size())
            return false;

        const_iterator seg1, seg2;

        for (seg1 = d1.begin(), seg2 = begin(); seg2 != end(); seg1++, seg2++) {
            if (seg1->first != seg2->first)
                return false;
            if (seg1->second != seg2->second)
                return false;
        }

        return true;
    }

    ReadManager::ReadManager() {
    }

    ReadManager::ReadManager(const string& filePath)
        : ifstream(filePath.c_str(), ios::in | ios::binary) {
        exceptions(ios::badbit | ios::eofbit | ios::failbit);
        // open(filePath);
        readHeader();
    }

    void ReadManager::open(const string& filePath) {
        exceptions(ios::badbit | ios::eofbit | ios::failbit);

        ifstream::open(filePath.c_str(), ios::in | ios::binary);
        readHeader();
    }

    void ReadManager::readHeader() {
        char readBuffer[HEADER_LENGTH + 1];

        seekg(ios::beg);

        ifstream::read(readBuffer, HEADER_LENGTH);

        readBuffer[HEADER_LENGTH] = 0;  // c string terminator
        fileInfo = readBuffer;

        cellInfoOrigin = tellg();

        isNative = fileInfo.sourceIsNative();

        CellInfo firstCell = retrieveCellInfo(0);

        startLocation = firstCell.getDataLocation();
        mappingLocation = firstCell.getMappingLocation();

        mappingSize_floats = fileInfo.getMappingSize() * fileInfo.getTotalNumberOfCompartments();
    }

    CellInfo ReadManager::makeNullInfo() {
        char fakeBuffer[SIZE_CELL_INFO_LENGTH];

        for (char *i = fakeBuffer, *bufferEnd = fakeBuffer + SIZE_CELL_INFO_LENGTH; i < bufferEnd;
             i++)
            *i = 0;

        return CellInfo(fakeBuffer);
    }

    const CellInfo ReadManager::_null_info = ReadManager::makeNullInfo();

    CellInfo ReadManager::retrieveFindCell(CellID cellNum) {
        int start = (cellNum > fileInfo.getNumberOfCells())
                        ? fileInfo.getNumberOfCells() - 1
                        : cellNum - 1;  // in many cases cellNum will be same as cellIndex
        CellInfo seek;

        bool oppositeSearch = false;

        int i = start;
        do {
            if (i >= (int)fileInfo.getNumberOfCells() || i < 0) {
                if (oppositeSearch)
                    return makeNullInfo();
                else {
                    oppositeSearch = true;
                    i = start;
                }
            }

            seek = retrieveCellInfo(i);

            if ((int)seek.getCellNum() < i)
                if (oppositeSearch)
                    i--;
                else
                    i++;
            else if (oppositeSearch)
                i++;
            else
                i--;

        } while (seek.getCellNum() != cellNum);

        return seek;
    }

    CellInfo ReadManager::retrieveCellInfo(CellIndex cellIndex, CellInfo::key sortMode) {
        FilePosition ppos = tellg();

        seekgCellInfo(cellIndex);

        char readBuffer[SIZE_CELL_INFO_LENGTH];
        read(readBuffer, SIZE_CELL_INFO_LENGTH);

        CellInfo reqInfo(readBuffer, isNative, sortMode);

        seekg(ppos);

        return reqInfo;
    }

    FrameInfo ReadManager::retrieveAllCellInfo(CellInfo::key sortMode) {
        FilePosition ppos = tellg();

        CellCount amountOfCells = fileInfo.getNumberOfCells();
        streamsize bytesToRead = SIZE_CELL_INFO_LENGTH * amountOfCells;
        char* readBuffer = new char[bytesToRead];
        char* endOfBuffer = readBuffer + bytesToRead;
        seekgCellInfo(0);

        read(readBuffer, bytesToRead);

        FrameInfo allCells(readBuffer, endOfBuffer, isNative, sortMode);

        delete[] readBuffer;

        seekg(ppos);

        return allCells;
    }

    void ReadManager::readFrame(DataItem* dataBuffer) {
        read(dataBuffer, fileInfo.getTotalNumberOfCompartments());
    }

    void ReadManager::readFrameMapping(MappingItem* buffer) {
        FilePosition ppos = tellg();

        seekgMappingStart();

        read(buffer, mappingSize_floats);

        seekg(ppos);
    }

    FrameParser::FrameOrganiser::FrameOrganiser() {
        sectSort = SID;
    }

    FrameParser::FrameOrganiser::FrameOrganiser(SortingKey segmentSorting) {
        sectSort = segmentSorting;
    }

    void FrameParser::FrameOrganiser::fillData(const FrameInfo& info,
                                               const vector<CompartmentMapping>& mapping,
                                               const DataItem* dataBuffer) {
        FrameInfo::const_iterator pinf = info.begin();
        vector<CompartmentMapping>::const_iterator pmap = mapping.begin();

        for (const DataItem *ptr = dataBuffer; pinf != info.end();
             ptr += pinf->getAmountOfCompartments(), pinf++, pmap++) {
            // iteration through all gids

            iterator ins = insert(make_pair(*pinf, CellData(sectSort)))
                               .first;         // insert CellData object for current gid
            ins->second.fillData(ptr, *pmap);  // feed new CellData object with mapping and data
        }
    }

    void FrameParser::FrameOrganiser::makeOrderedDataBuffer(DataItem* buffer) const {
        for (const_iterator cd = begin(); cd != end(); cd++)  // iterate through cells
            for (CellData::const_iterator sd = cd->second.begin(); sd != cd->second.end();
                 sd++)  // iterate through segments
                for (SectionData::const_iterator f = sd->second.begin(); f != sd->second.end();
                     f++, buffer++)  // iterate through compartments
                    *buffer = *f;
    }

    FrameParser::CompartmentIndexing FrameParser::FrameOrganiser::makeOrderedOffsetReferences()
        const {
        CompartmentIndexing gidV;
        gidV.reserve(size());

        FrameDataIndex position = 0;

        for (const_iterator gid = begin(); gid != end(); gid++) {  // iterate through gids
            gidV.push_back(vector<FrameDataIndex>());

            // reserve memory for amount of segments
            gidV.back().reserve((int)(gid->second.size() * 1.5));

            MappingItem presid = 0;
            for (CellData::const_iterator sid = gid->second.begin(); sid != gid->second.end();
                 position += sid->second.size(), sid++) {
                // iterating through available sid's
                while (presid++ != sid->first.getSID())
                    gidV.back().push_back(0);  // insert dummy variables for unavailable sids (gaps)
                gidV.back().push_back(position);
            }
        }

        return gidV;
    }

    FrameParser::CompartmentCounts FrameParser::FrameOrganiser::makeOrderedSegSizes() const {
        CompartmentCounts gidV;
        gidV.reserve(size());

        for (const_iterator gid = begin(); gid != end(); gid++) {
            gidV.push_back(vector<SectionDataCount>());
            gidV.back().reserve((int)(gid->second.size() * 1.5));
            MappingItem presid = 0;
            for (CellData::const_iterator sid = gid->second.begin(); sid != gid->second.end();
                 sid++) {
                while (presid++ != sid->first.getSID())
                    gidV.back().push_back(0);  // fill sid gaps with zero sizes
                gidV.back().push_back(sid->second.size());
            }
        }

        return gidV;
    }

    FrameParser::FrameParser() {
        readBuffer = NULL;
        refArray = NULL;
    }

    FrameParser::FrameParser(const string& file, bool sortData) : ReadManager(file) {
        if (fileInfo.getMappingSize() != 1) {
            cerr
                << "Target file has incompatible mapping size for FrameParser. FrameParser only works with mapping size 1."
                << endl;
            exit(1);
        }

        CellInfo::key cinfoSorting = (sortData) ? CellInfo::GID : CellInfo::LOCATION;

        FrameInfo cinfo = retrieveAllCellInfo(cinfoSorting);

        ReadOptimiser::setDefaultElementSize(sizeof(DataItem));
        offsets = ReadOptimiser(fileInfo.getTotalNumberOfCompartments());  // * sizeof(DataItem) );
        createReferences(cinfo, sortData);

        timestep = 0;  // by default start pointing to timstep 0
        ReadManager::seekgFrame(timestep);
    }

    FrameParser::FrameParser(const string& file, const vector<CellID>& gidTarget, bool sortData)
        : ReadManager(file) {
        timestep = 0;
        ReadOptimiser::setDefaultElementSize(sizeof(DataItem));
        retarget(gidTarget, sortData);
    }

    FrameParser::~FrameParser() {
        if (readBuffer)  // readBuffer may have not been used if sorting was disabled, in which case
                         // it would have been set to NULL
            delete[] readBuffer;
        if (refArray)  // same applies for refArray
            delete[] refArray;
    }

    void FrameParser::retarget(const vector<CellID>& gidTarget, bool sortData) {
        CellInfo::key cinfoSorting = (sortData) ? CellInfo::GID : CellInfo::LOCATION;

        FrameInfo cinfo = retrieveAllCellInfo(cinfoSorting);
        cinfo.filterGIDs(gidTarget);

        offsets = ReadOptimiser(cinfo, getFrameSize_bytes(), getStartOffset());
        offsets.setElementSize(sizeof(DataItem));

        createReferences(cinfo, sortData);

        seekg(getStartOffset() + offsets.getFirst());

        operator=(timestep);
    }

    void FrameParser::createReferences(const FrameInfo& cinfo, bool sortData) {
        SortingKey sSorting = (sortData) ? SID : LOC;

        refArray = NULL;
        readBuffer = NULL;

        vector<CompartmentMapping> allMapping;

        DataItem* tempBuffer =
            new DataItem[offsets.getDataSize_elements()];  // allocate memory for readBuffer. This
                                                           // will atleast be needed for creating
                                                           // the FrameOrganiser

        readFrameMapping(tempBuffer);

        FrameInfo::const_iterator pinf = cinfo.begin();

        for (DataItem *mpptr = tempBuffer, *endptr = mpptr + pinf->getAmountOfCompartments();
             pinf != cinfo.end(); pinf++, mpptr = endptr, endptr += pinf->getAmountOfCompartments())
            allMapping.push_back(CompartmentMapping(mpptr, endptr));

        if (sortData)
            for (FrameDataIndex i = 0; i < offsets.getDataSize_elements(); i++)
                tempBuffer[i] = i;  // makes buffer of indexes { 0 , 1 ,2 , 3 , ... }
        // else it doesnt really matter what FrameOrganiser is fed with.

        org = FrameOrganiser(sSorting);

        org.fillData(cinfo, allMapping, tempBuffer);  // feed the FrameOrganiser with the data.

        // references = org.makeOrderedOffsetReferences();
        // segSizes = org.makeOrderedSegSizes();

        if (sortData) {
            readBuffer = tempBuffer;  // keep allocated memory at tempBuffer

            // sortData is enabled -> create an array of pointers to readBuffer
            // in which the pointers are sorted by GID and SID of the data the point to.
            // these pointers are to be stored into ref array
            org.makeOrderedDataBuffer(
                readBuffer);  // rearranges indexes to correspond to gid and segment order
            // the memory allocated for readBuffer is now being used for a different purpose:
            // it contains the sorted indexes, which now need to be converted to pointers to
            // readBuffer.

            refArray = new DataItem*[offsets.getDataSize_elements()];  // allocate memory for
                                                                       // refArray (will contain
                                                                       // pointers for sorting)

            for (float *pIndex = readBuffer, **pPtr = refArray,
                       *pEnd = readBuffer + offsets.getDataSize_elements();
                 pIndex != pEnd; pIndex++, pPtr++)
                *pPtr = &readBuffer[(FrameDataIndex)*pIndex];  // converts indexes to pointers and
                                                               // stores them into refArray
        } else
            // data sorting is disabled -> allocated memory at temp Buffer is no longer needed.
            //(data will be read directly into the memory that user has provided)
            delete[] tempBuffer;
    }

    void FrameParser::bufferTransfer(DataItem* buffer) const {
        DataItem* reqLocation = buffer;

        if (buffer == readBuffer)  // source and destination is the same - trouble!
            buffer = new DataItem[offsets.getDataSize_elements()];  // need to allocate different
                                                                    // destination

        for (DataItem **sourcePtr = refArray, *target = buffer,
                      *dEnd = buffer + offsets.getDataSize_elements();
             target != dEnd; sourcePtr++, target++)
            *target = **sourcePtr;

        if (reqLocation != buffer) {  // buffer has not been transfered where requested!
            // this can now be corrected.
            for (DataItem *origin = buffer, *destination = reqLocation,
                          *dEnd = buffer + offsets.getDataSize_elements();
                 origin != dEnd; origin++, destination++)
                *destination = *origin;  // copy data to requested location
            delete[] buffer;  // buffer was is only a temporary memory allocation - delete it.
        }
    }

    FrameParser& FrameParser::operator++() {
        timestep++;
        seekg(fileInfo.getTotalNumberOfCompartments() * sizeof(DataItem), cur);
        return *this;
    }

    FrameParser& FrameParser::operator++(int) {
        return this->operator++();
    }

    FrameParser& FrameParser::operator--() {
        timestep--;
        seekg(-fileInfo.getTotalNumberOfCompartments() * sizeof(DataItem), cur);
        return *this;
    }

    FrameParser& FrameParser::operator--(int) {
        return this->operator--();
    }

    FrameParser& FrameParser::operator+=(FrameDiff increment) {
        timestep += increment;
        seekg(increment * fileInfo.getTotalNumberOfCompartments() * sizeof(DataItem), cur);
        return *this;
    }

    FrameParser& FrameParser::operator-=(FrameDiff decrement) {
        timestep -= decrement;
        seekg(-decrement * fileInfo.getTotalNumberOfCompartments() * sizeof(DataItem), cur);
        return *this;
    }

    FrameParser& FrameParser::operator=(FrameIndex newTime) {
        timestep = newTime;
        seekgFrame(timestep);
        seekg(offsets.getFirst(), cur);
        return *this;
    }

    Time FrameParser::getTime() const {
        return fileInfo.getStartTime() + timestep * fileInfo.getTimeStepSize();
    }

    void FrameParser::readFrame(DataItem* buffer) {
        DataItem* readPtr;

        if (readBuffer)            // readBuffer will be NULL if sorting is disabled
            readPtr = readBuffer;  // sorting enabled -> read into readBuffer and transfer data in
                                   // order later
        else
            readPtr = buffer;  // sorting disabled -> read directly into buffer

        // offsets contains alternating values of data sizes to be read and skipped.
        for (ReadOptimiser::iterator i = offsets.begin(); i != offsets.end(); i++) {
            read(readPtr, *i);
            readPtr += *i;
            i++;

            // seems to be confusion as to whether offsets should contain bytelength or element
            // count.  I have it contain element count for now
            // i now points to an offset that indicates amount of data to be skipped
            if (i == offsets.end())  // no more offsets -> i'm done reading
                break;

            if (*i)
                seekg(*i, ios::cur);  // non-zero offset -> perfom the skip
        }

        if (refArray)                // refArray will be NULL if sorting is disabled
            bufferTransfer(buffer);  // sorting is enabled. Transfer data from readBuffer to buffer
                                     // such that data is in order
    }

    void FrameParser::readFrameMapping(MappingItem* buffer) {
        // readFrameMapping should not return without the filepointer being in the same state as
        // before it was called.
        FilePosition before = tellg();  // make a backup for current filepointer

        seekgMappingStart();  // now move file pointer to mapping location

        size_t prevElement = offsets.getElementSize();  // backup element size.
        offsets.setElementSize(sizeof(MappingItem));    // setup offsets to be reading MappingItems
        // instead of whatever it was reading before.

        seekg(offsets.getFirst(), ios::cur);  // and advance to first cell
        readFrame(
            buffer);  // read the data - this function will also alter the filepointer and timestep

        offsets.setElementSize(prevElement);  // restore element size.

        seekg(before);  // restore filepointer
    }

    void FrameParser::readFrameData(DataItem* buffer) {
        readFrame(buffer);
        timestep++;  // the filepointer has been incremented to the next timestep as data has been
                     // read.
    }

    size_t FrameParser::ReadOptimiser::defaultElementSize = 1;

    FrameParser::ReadOptimiser::ReadOptimiser() {
        elementSize = 0;
        firstOffset = 0;
        dataSize_elements = 0;
    }

    FrameParser::ReadOptimiser::ReadOptimiser(std::streamsize frameSize_bytes) {
        elementSize = defaultElementSize;
        firstOffset = 0;
        dataSize_elements = frameSize_bytes;  // * elementSize;
        push_back(dataSize_elements);
    }

    FrameParser::ReadOptimiser::ReadOptimiser(const FrameInfo& cellsToRead,
                                              std::streamsize frameSize_bytes,
                                              FilePosition startLocation) {
        elementSize = defaultElementSize;
        // collect data edge points
        multiset<FilePosition> dataEdges;

        dataEdges.insert(startLocation);

        for (FrameInfo::const_iterator i = cellsToRead.begin(); i != cellsToRead.end(); i++) {
            dataEdges.insert(i->getDataLocation());  // starting point
            dataEdges.insert(i->getDataLocation() +
                             i->getAmountOfCompartments() * elementSize);  // end point
        }

        dataEdges.insert(startLocation + (FileOffset)frameSize_bytes);

        /*create a vector that contains the offsets.
         * The first offset represents the first bunch of data to be read.
         * The second offset then represents the amound of data to be skipped after the first read.
         * and then the next offset is again the amount of data to be read.
         */

        // compute the offsets differentiating the dataEdges
        list<FileOffset> offsets;
        for (multiset<FilePosition>::iterator currentPos = dataEdges.begin(),
                                              prevPos = currentPos++;
             currentPos != dataEdges.end(); prevPos = currentPos, currentPos++)
            offsets.push_back(*currentPos - *prevPos);

        // the last offsets should move the file pointer beyond the frame to the start location on
        // the next frame
        // therefore the first and last offsets are special -> take them out

        firstOffset = offsets.front();
        FileOffset lastOffset = offsets.back();

        offsets.pop_front();  // remove first offset from list
        offsets.pop_back();   // remove last offset from list

        // now the first an last offset are read dimensions. This guarantees that a skip has always
        // a precessor and successor
        // this is convenient in the next step, as a skip is always on an even location aswell as
        // offsets.end()

        list<FileOffset>::iterator it;
        bool toggleOnRead;

        dataSize_elements = 0;
        for (toggleOnRead = true, it = offsets.begin(); it != offsets.end();
             toggleOnRead = !toggleOnRead, it++) {
            if (toggleOnRead) {
                dataSize_elements += (*it /= elementSize);
                // handle read
                // read specifications should be in units of floats, but skips should be in units of
                // bytes.
                // also size of data to be read per frame needs to established
            } else
                // handle skip - There should be no zero offset inbetween
                if (*it == 0) {
                // remove this zero skip by combining prev. and next read.
                list<FileOffset>::iterator prevRead = it;
                list<FileOffset>::iterator nextRead = it;
                prevRead--;
                nextRead++;

                dataSize_elements +=
                    (*nextRead /= elementSize);  // next read was still in units of bytes
                *prevRead += *nextRead;          // combine nextRead with prevRead
                offsets.erase(it);
                offsets.erase(nextRead);
                it = prevRead;  // point to prevRead now
                toggleOnRead =
                    true;  // and mark that toggle flag that current iteration deal with read.
            }
        }

        // now optimise access performance by converting the offsets list to a vector.

        vector<FileOffset> voffsets(offsets.begin(), offsets.end());
        swap(voffsets);

        // compute the finalSkip, which brings the filePointer back to the starting location, but of
        // the nextFrame
        FileOffset finalSkip = lastOffset + firstOffset;

        if (finalSkip)
            push_back(finalSkip);
    }

    void FrameParser::ReadOptimiser::updateOffsets(size_t newSize) {
        dataSize_elements = newSize * dataSize_elements / elementSize;
        for (iterator i = begin(); i != end(); i++) {  // iterate through all read offsets
            // first reduce the read offsets to byte numbers, then convert it to the new element
            // size
            *i = newSize * (*i) / elementSize;
            i++;
            if (i == end())
                break;
        }
    }
}
