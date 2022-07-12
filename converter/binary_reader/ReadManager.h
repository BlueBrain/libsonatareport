#ifndef _ReadManager__
#define _ReadManager__

#include <stdint.h>

#include <fstream>
#include <string>

#include <vector>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <stdexcept>

#include "BinaryPositions.h"

#include "crossPlatformConversion.h"

/*Terminology:
 * Compartment: lowest subdivision of cell. (cannot be broken down further)
 * Section: group of unbranched compartments (may sometimes be called segment, please change to
 * section if found)
 * Frame: Data residing in a single timestep. The scope of a frame may vary in terms of which cells
 * it includes, depending on its use.
 * See documentation for further naming conventions.
 */

namespace bbpReader {

    // class list:
    class Header;
    class CellInfo;
    class FrameInfo;
    class CompartmentMapping;
    class SectionData;
    // class SectionData::Identifier;
    class Identifier;
    class CellData;
    class ReadManager;
    class FrameParser;
    // class FrameParser::FrameOrganiser;
    // class FrameParser::ReadOptimiser;
    // end class list

    // type conventions:

    typedef uint32_t CellCount;   // indications amount of cells.
    typedef uint32_t CellID;      // used for GIDs
    typedef uint32_t CellIndex;   // used for indexing arrays of cells.s
    typedef uint32_t FrameCount;  // indictions of amount of timesteps.
    typedef uint32_t FrameIndex;  // used for frame numbers
    typedef int32_t
        FrameDiff;  // used for indicating time periods in terms of frame indexes (can be negative)
    typedef double Time;              // used for indicating biological time durations.
    typedef uint64_t FrameDataCount;  // used for indicating amount of data items for an entire
                                      // frame
    typedef uint64_t FrameDataIndex;  // used for indexing arrays, vectors etc of data items.
    typedef uint32_t CellDataCount;   // used for indicating amount of data items that reside inside
                                      // a single Cell.
    typedef uint16_t SectionDataCount;  // used for indicating amount of data items that reside
                                        // inside a single section.
    typedef uint16_t SectionIndex;      // used for indexing data inside single sections.
    typedef float MappingItem;  // type in which SID and other mapping types are represented in.
    typedef int32_t ExtraMappingItem;  // type in which extra mapping items are represented in.
    typedef float DataItem;  // type in which the actual voltage, current, conductivity etc values
                             // are represented in

    enum SortingKey { SID, LOC };  // sorting specifier for sections.
    // SID will make sections to be sorted by their mapping,
    // LOC by their locaiton in file.

    namespace baseTypes {
        // these typedefs define from which STL some of the clases are derived from

        typedef std::set<bbpReader::CellInfo> FrameInfo;
        typedef std::vector<MappingItem> CompartmentMapping;
        typedef std::vector<DataItem> SectionData;
        // typedef std::map<bbpReader::SegmentData::Identifier,bbpReader::SegmentData>
        // CellData; done blow
        typedef std::map<bbpReader::CellInfo, bbpReader::CellData> FrameOrganiser;
        typedef std::vector<std::ios::off_type> ReadOptimiser;
    }

    // end of type conventions

    class Header {
      private:
        bool isNative;

        // header info
        uint32_t headerSize;
        std::string libraryVersion;
        std::string simulatorVersion;
        CellCount amountOfCells;
        FrameDataCount totalAmountOfCompartments;
        FrameCount numberOfSteps;
        Time startTime;
        Time endTime;
        Time dt;
        std::string dataUnit;
        std::string timeUnit;
        FrameDataCount mappingSize;
        std::string mappingName;
        FrameDataCount extraMappingSize;
        std::string extraMappingName;
        std::string targetName;

      public:
        Header() : headerSize(0){}
        Header(char* readBuffer);
        // Header(const Header& source);

        bool operator==(const Header& h1) const;

        // get methods:
        inline bool sourceIsNative() const {
            return isNative;
        }
        inline std::string getLibraryVersion() const {
            return libraryVersion;
        }
        inline std::string getSimulatorVersion() const {
            return simulatorVersion;
        }
        inline CellCount getNumberOfCells() const {
            return amountOfCells;
        }
        inline FrameDataCount getTotalNumberOfCompartments() const {
            return totalAmountOfCompartments;
        }
        inline FrameCount getNumberOfSteps() const {
            return numberOfSteps;
        }
        inline std::string getMappingName() const {
            return mappingName;
        }
        inline std::string getExtraMappingName() const {
            return extraMappingName;
        }
        inline std::string getDataUnit() const {
            return dataUnit;
        }
        inline std::string getTimeUnit() const {
            return timeUnit;
        }
        inline FrameDataCount getMappingSize() const {
            return mappingSize;
        }
        inline FrameDataCount getExtraMappingSize() const {
            return extraMappingSize;
        }
        inline std::string getTargetName() const {
            return targetName;
        }
        inline Time getStartTime() const {
            return startTime;
        }
        inline Time getEndTime() const {
            return endTime;
        }
        inline Time getTimeStepSize() const {
            return dt;
        }
        // end get methods
    };

    class CellInfo {
      private:
        CellID cellNum;
        CellDataCount amountOfCompartments;

        std::ios::off_type dataLocation;
        std::ios::off_type extraMappingLocation;
        std::ios::off_type mappingLocation;

        void* master;

      public:
        enum key { GID, LOCATION };

        CellInfo();
        CellInfo(char* readBuffer, bool isNative = true, key cmp = GID);

        bool operator==(const CellInfo& c1) const;
        bool operator!=(const CellInfo& c1) const;
        bool operator<=(const CellInfo& c1) const;
        bool operator>=(const CellInfo& c1) const;
        bool operator<(const CellInfo& c1) const;
        bool operator>(const CellInfo& c1) const;

        inline std::ios::off_type getDataLocation() const {
            return dataLocation;
        }

        inline std::ios::off_type getMappingLocation() const {
            return mappingLocation;
        }

        inline std::ios::off_type getExtraMappingLocation() const {
            return extraMappingLocation;
        }

        inline CellID getCellNum() const {
            return cellNum;
        }

        inline CellDataCount getAmountOfCompartments() const {
            return amountOfCompartments;
        }

        static bool compareByOffset(const CellInfo& c1, const CellInfo& c2);
    };

    class FrameInfo : private baseTypes::FrameInfo {
      private:
        typedef baseTypes::FrameInfo Base;
        typedef std::map<CellID, iterator> Lookup;
        Lookup gidLookup;

      public:
        typedef Base::const_iterator const_iterator;

        FrameInfo() {
        }
        FrameInfo(char* buffer,
                  char* endOfBuffer,
                  bool isNative = true,
                  CellInfo::key sortMode = CellInfo::GID);
        FrameInfo& filterGIDs(const std::vector<CellID>& gids);

        /**
         * Copy constructor.  Needed to ensure set iterators reference the correct collection
         */
        FrameInfo(const FrameInfo& rhs);

        /**
         * Assignment operator.  Needed to ensure set iterators reference the correct collection
         */
        FrameInfo& operator=(const FrameInfo& rhs);

        inline const CellInfo& operator[](CellID gid) {
            return *gidLookup[gid];
        }

        inline const_iterator begin() const {
            return Base::begin();
        }

        inline const_iterator end() const {
            return Base::end();
        }
    };

    class CompartmentMapping : public baseTypes::CompartmentMapping {
      private:
        typedef baseTypes::CompartmentMapping Base;

      public:
        typedef Base::size_type SizeType;

        CompartmentMapping();
        CompartmentMapping(MappingItem* start, MappingItem* end);
        CompartmentMapping(MappingItem* buffer, SizeType bufSize);
    };

    class SectionData : public baseTypes::SectionData {
      private:
        typedef baseTypes::SectionData Base;

      public:
        class Identifier;
        SectionData();
        SectionData(DataItem firstData);
        SectionData(DataItem* start, DataItem* end);

        bool operator==(const SectionData& s1);
        bool operator!=(const SectionData& s1);

        inline SectionDataCount numOfCompartments() {
            return size();
        }

        static bool compareByOriginalLocation(const Identifier& c1, const Identifier& c2);
    };

    class SectionData::Identifier {
      public:
      private:
        typedef SectionData::Identifier thisType;

        MappingItem sid;
        SectionIndex loc;

        SortingKey key;

      public:
        inline Identifier(const MappingItem& map,
                          const SectionIndex& location,
                          SortingKey primary = SID) {
            sid = map, loc = location, key = primary;
        }

        inline bool operator!=(const thisType& c) const {
            return (key == SID) ? sid != c.sid : loc != c.loc;
        }

        inline bool operator==(const thisType& c) const {
            return (key == SID) ? sid == c.sid : loc == c.loc;
        }

        inline bool operator<(const thisType& c) const {
            return (key == SID) ? sid < c.sid : loc < c.loc;
        }

        inline bool operator<=(const thisType& c) const {
            return (key == SID) ? sid <= c.sid : loc <= c.loc;
        }

        inline bool operator>(const thisType& c) const {
            return (key == SID) ? sid > c.sid : loc > c.loc;
        }

        inline bool operator>=(const thisType& c) const {
            return (key == SID) ? sid >= c.sid : loc >= c.loc;
        }

        inline operator MappingItem() const {
            return sid;
        }

        inline MappingItem getSID() const {
            return sid;
        }

        inline SectionIndex getLocation() const {
            return loc;
        }
    };

    namespace baseTypes {
        typedef std::map<bbpReader::SectionData::Identifier, bbpReader::SectionData> CellData;
    }

    class CellData : public baseTypes::CellData {
      private:
        typedef baseTypes::CellData Base;
        typedef std::map<MappingItem, iterator> Lookup;

        SortingKey sectSort;
        Lookup sidLookup;

      public:
        CellData();
        CellData(const SortingKey segmentSorting = SID);
        CellData(const DataItem* dataBuffer, const CompartmentMapping& mapping);
        void fillData(const DataItem* dataBuffer, const CompartmentMapping& mapping);

        inline iterator find(const SectionData::Identifier& i) {
            return Base::find(i);
        }

        inline iterator find(const MappingItem& map_find) {
            Lookup::iterator match = sidLookup.find(map_find);
            if (match == sidLookup.end())
                return end();
            else
                return match->second;
        }

        bool operator==(const CellData& d1);
    };

    class ReadManager : public std::ifstream {
      private:
        typedef std::ifstream Base;

      public:
        typedef Base::off_type FileOffset;
        typedef Base::pos_type FilePosition;

      private:
        FileOffset cellInfoOrigin;  // position in file where cell info starts (after header)

        bool isNative;  // indicates whether file is originated from an architecture that arranges
                        // variables in same byte order

        FrameDataCount mappingSize_floats;

        FileOffset startLocation;
        FileOffset mappingLocation;

        inline void seekgStart() {
            seekg(startLocation, beg);
        }

        static CellInfo makeNullInfo();

      protected:
        inline FileOffset getStartOffset() {
            return startLocation;
        }

      public:
        static const CellInfo _null_info;

        Header fileInfo;

        ReadManager();

        /**
         * Constructor which immediately opens a file for reading.
         *
         * @param filePath The bbp binary report file to open
         */
        explicit ReadManager(const std::string& filePath);

        /**
         * Opens the filename specified provided this object has not already opened a file
         *
         * @param filePath The bbp binary report file to open
         */
        void open(const std::string& filePath);

        // Find the cell's info which corresponds to cellNum (note: cellNum is not the cell index in
        // the file)
        CellInfo retrieveFindCell(CellID cellNum);

        // Retrieve  global (time invariant) information of a cell (not actual data of the cell)
        CellInfo retrieveCellInfo(CellIndex cellIndex, CellInfo::key sortMode = CellInfo::LOCATION);

        // NOTE: its highly recommended to either disable cache (set to 0) or set cache size to
        // unlimited (set to -1) before using this method:
        FrameInfo retrieveAllCellInfo(CellInfo::key sortMode = CellInfo::LOCATION);

        /*seekg*:
         * all the seekg functions below are wrapping around the original iostream seekg method.
         */

        inline void seekgCellInfo(CellIndex cellIndex = 0) {
            seekg(cellInfoOrigin);
            if (cellIndex)
                seekg(cellIndex * SIZE_CELL_INFO_LENGTH, cur);
        }

        inline void seekgData(const CellInfo& cellSpec, FrameIndex timeStep = 0) {
            seekg(cellSpec.getDataLocation() + timeStep * getFrameSize_bytes(), beg);
        }

        inline void seekgMapping(const CellInfo& cellSpec) {
            seekg(cellSpec.getMappingLocation(), beg);
        }

        inline void seekgExtraMapping(const CellInfo& cellSpec) {
            seekg(cellSpec.getExtraMappingLocation(), beg);
        }

        inline void seekgMappingStart() {
            seekg(mappingLocation, beg);
        }

        inline void seekgFrame(FrameIndex timeStep) {
            seekg(startLocation + getFrameSize_bytes() * timeStep);
        }

        /*read():
         * this function essentially takes 3 arguments (the third being the datatype).
         * it does the same as the read function of the base ifstream class that it is overloading
         * the only difference is that num is not neccesarily in number of bytes
         * instead it defines the mount of elemetns in the buffer (depending on the datatype of the
         * buffer).
         *
         * Important:
         * This read method will automatically fix the byte orders of the data if applicable.
         * The algorithm of fixing the byte orders depends on the size of each element in the buffer
         * array.
         * Therefore you should be carefull with passing char* buffers,
         * unless you are really intending to read single byte data.
         */

        template <class datatype>
        ReadManager& read(datatype* buffer, std::streamsize num) {
            Base::read((char*)buffer, num * sizeof(datatype));

            if (sizeof(datatype) != 1 && !isNative)  // a good compiler should not include this
                                                     // entire if statement when this function is
                                                     // reading single byte items.
                // This if is to avoid entering an empty loop. Otherwise could have been done in
                // switch below.
                for (datatype* endOfBuffer = buffer + num; buffer < endOfBuffer; buffer++)
                    switch (sizeof(datatype)) {  // metamorphic switch
                        case 8:                  // double, int64_t, offsets etc
                            DOUBLESWAP(*buffer);
                            break;
                        case 4:  // float, int32_t, etc
                            FIX_INT(*buffer);
                            break;
                        case 2:  // short
                            FIX_SHORT(*buffer);
                            break;
                        default:  // anything weird that cannot be fixed.
                            // TODO: throw error exception
                            break;
                    }
            return *this;
        }

        /*checkState():
         * will determine whether a reading error has occured,
         * output the according error message if there has been an error,
         * and exit execution if an error occured.
         */

        void checkState();

        /*getFrameSize:
         * returns the number of bytes that hold the data of all compartments of all cells for a
         * single timestep.
         * please note that its is in bytes, not floats!
         * use fileInfo.getTotalNumberOfCompartments() instead for number of floats.
         */

        inline std::streamsize getFrameSize_bytes() const {
            return fileInfo.getTotalNumberOfCompartments() * sizeof(DataItem);
        }

        /*getMappingSize:
         * returns the number of bytes that hold the data of the mappings of all compartmetns of all
         * cells.
         * please note that this is in bytes not floats.
         */

        inline std::streamsize getMappingSize_bytes() const {
            return mappingSize_floats * sizeof(MappingItem);
        }

        // The following are low level c-style reading interfaces for max performance - they provide
        // unsorted, raw data.
        /*readFrameMapping:
         * will fill buffer with mappings of all compartmetns of all cells.
         * the file pointer will be automatically adjusted accordingly and then restored to original
         * position.
         * use getMappingSize() / sizeof(float) to determine the size of the buffer in number bytes
         * (NOT float).
         */

        void readFrameMapping(MappingItem* buffer);

        /*readFrame:
         * will fill dataBuffer with data from next frame.
         * use seekgFrame() to specify timestep.
         * use fileInfo.getTotalNumberOfCompartments() to determine the size of buffer
         */
        void readFrame(DataItem* dataBuffer);

      private:
        /**
         * Upon opening a report file, read the header with metadata
         */
        void readHeader();
    };

    class FrameParser : protected ReadManager {
      public:
        typedef std::vector<std::vector<FrameDataIndex> > CompartmentIndexing;
        typedef std::vector<std::vector<SectionDataCount> > CompartmentCounts;

      protected:
        class FrameOrganiser : public baseTypes::FrameOrganiser {
          private:
            SortingKey sectSort;

          public:
            FrameOrganiser();
            explicit FrameOrganiser(SortingKey sectionSorting);

            void fillData(const FrameInfo& info,
                          const std::vector<CompartmentMapping>& mapping,
                          const float* dataBuffer);

            void makeOrderedDataBuffer(DataItem* buffer) const;

            CompartmentIndexing makeOrderedOffsetReferences() const;
            CompartmentCounts makeOrderedSegSizes() const;
        };

        class ReadOptimiser : public baseTypes::ReadOptimiser {
          private:
            void updateOffsets(size_t newSize);

          protected:
            FileOffset firstOffset;
            FrameDataCount dataSize_elements;

            size_t elementSize;
            static size_t defaultElementSize;

          public:
            ReadOptimiser();
            explicit ReadOptimiser(std::streamsize frameSize_bytes);
            ReadOptimiser(const FrameInfo& cellsToRead,
                          std::streamsize frameSize_bytes,
                          FilePosition startLocation);

            inline FileOffset getFirst() const {
                return firstOffset;
            }

            inline FrameDataCount getDataSize_elements() const {
                return dataSize_elements;
            }

            inline void setElementSize(size_t newSize) {
                if (newSize != elementSize) {
                    updateOffsets(newSize);
                    elementSize = newSize;
                }
            }

            static inline void setDefaultElementSize(size_t newDefault) {
                defaultElementSize = newDefault;
            }

            inline size_t getElementSize() const {
                return elementSize;
            }
        };

        void readFrame(DataItem* buffer);

      private:
        DataItem* readBuffer;
        DataItem** refArray;

        ReadOptimiser offsets;
        FrameOrganiser org;

        FrameIndex timestep;

        void createReferences(const FrameInfo& cinfo, bool sortData);

        void bufferTransfer(DataItem* target) const;

      public:
        /*Constructor:
         * Will read next frame of data and parse it into a map of CellData structure from file.
         */
        FrameParser();
        FrameParser(const std::string& file, bool sortData = true);
        FrameParser(const std::string& file,
                    const std::vector<CellID>& gidTarget,
                    bool sortData = true);
        ~FrameParser();

        /*retarget:
         * Takes a list of gids in form of a vector of CellID and uses it to respecify the target
         * cells of this
         * frameparser.
         */

        inline FrameIndex simtime2index(Time time) {
            return (FrameIndex)((time - fileInfo.getStartTime()) / fileInfo.getTimeStepSize());
        }

        void retarget(const std::vector<CellID>& gidTarget, bool sortData = true);

        inline const Header& getHeader() const {
            return fileInfo;
        }

        inline FrameDataCount getBufferSize_elements() const {
            return offsets.getDataSize_elements();
        }

        // to select frame (similar use to iterators):
        // Note that readFrameData already incorporates operator++ !
        FrameParser& operator++();                     // increment timestep
        FrameParser& operator++(int);                  // increment timestep
        FrameParser& operator--();                     // decrement timestep
        FrameParser& operator--(int);                  // decrement timestep
        FrameParser& operator+=(FrameDiff increment);  // increment timestep by scalar
        FrameParser& operator-=(FrameDiff decrement);  // decrement timestep by scalar
        FrameParser& operator=(FrameIndex newTime);    // set timestep to scalar

        // returns the current timestep:
        inline FrameIndex getTimestep() const {
            return timestep;
        }

        inline bool hasMore() const {
            return timestep < fileInfo.getNumberOfSteps();
        }

        // returns the current time that the current timestep represents.
        Time getTime() const;

        void readFrameMapping(MappingItem* buffer);
        void readFrameData(DataItem* buffer);  // reads data of next frame and advances to the frame
                                               // after. Data is sorted by GID and segment num

        /*getReferences()
         * returns a vector of vectors of offsets that correspond to the relevant segments inside
         * the read buffer.
         * eg.
         * getReferences()[12][31]; //will give the offset for segment 31 in cell 12
         * if the segment has no comparments, it will return an offset of 0.
         */
        inline CompartmentIndexing getReferences() const {
            return org.makeOrderedOffsetReferences();
        }

        /*getCompartmentCounts()
         * returns a vector of vectors like getReferences() above,
         * but instead of containing the offsets for each segment,
         * it contains the amount of compartmetns of the corresponding segment.
         */
        inline CompartmentCounts getCompartmentCounts() const {
            return org.makeOrderedSegSizes();
        }
    };
}
#endif
