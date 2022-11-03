#include "binary_reader.h"
#include <algorithm>
#include <iostream>

using namespace std;

namespace bbp {
namespace binary_reader {

Header::Header(char* read_buffer) {
#define fetch(TYPE, POSITION) *((TYPE*) (read_buffer + POSITION))
#define fetchS(POSITION) string((char*) (read_buffer + POSITION))

    double read_identifier = fetch(double, SONATA_REPORT_IDENTIFIER_POSITION);
    header_size = fetch(int32_t, SONATA_REPORT_HEADER_SIZE_POSITION);
    amount_of_cells = fetch(int32_t, SONATA_REPORT_TOTAL_NUMBER_OF_CELLS_POSITION);
    /** \bug 32 bit is consistent with the binary format specification of the
        header, however
     */
    uint32_t compartments_per_frame = fetch(int32_t,
                                            SONATA_REPORT_TOTAL_NUMBER_OF_COMPARTMENTS_POSITION);
    library_version = fetchS(SONATA_REPORT_LIBRARY_VERSION_POSITION);
    simulator_version = fetchS(SONATA_REPORT_SIMULATOR_VERSION_POSITION);
    number_of_steps = fetch(int32_t, SONATA_REPORT_NUMBER_OF_STEPS_POSITION);
    start_time = fetch(double, SONATA_REPORT_TIME_START_POSITION);
    end_time = fetch(double, SONATA_REPORT_TIME_END_POSITION);
    dt = fetch(double, SONATA_REPORT_DT_TIME_POSITION);
    data_unit = fetchS(SONATA_REPORT_D_UNIT_POSITION);
    time_unit = fetchS(SONATA_REPORT_T_UNIT_POSITION);
    mapping_size = fetch(int32_t, SONATA_REPORT_MAPPING_SIZE_POSITION);
    mapping_name = fetchS(SONATA_REPORT_MAPPING_NAME_POSITION);
    extra_mapping_size = fetch(int32_t, SONATA_REPORT_EXTRA_MAPPING_SIZE_POSITION);
    extra_mapping_name = fetchS(SONATA_REPORT_EXTRA_MAPPING_NAME_POSITION);
    target_name = fetchS(SONATA_REPORT_TARGET_NAME_POSITION);

#undef fetchS
#undef fetch

    is_native = read_identifier == SONATA_REPORT_ARCHITECTURE_IDENTIFIER;

    if (!is_native) {
        // check if data can be converted:
        DOUBLESWAP(read_identifier);
        if (read_identifier != SONATA_REPORT_ARCHITECTURE_IDENTIFIER) {
            cerr << "File is corrupt or originated from an unknown architecture." << endl;
            exit(1);
        }

        // if data was stored on a machine with different architecture,
        // then the data needs to be converted to be compatible with the reading machine
        FIX_INT(header_size);
        FIX_INT(amount_of_cells);
        FIX_INT(compartments_per_frame);
        FIX_INT(number_of_steps);
        DOUBLESWAP(start_time);
        DOUBLESWAP(end_time);
        DOUBLESWAP(dt);
        FIX_INT(mapping_size);
        FIX_INT(extra_mapping_size);
    }

    total_amount_of_compartments = compartments_per_frame;
}

bool Header::operator==(const Header& h1) const {
#define _CHECK(field)      \
    if (h1.field != field) \
    return false
    // compare all fields other than is_native
    _CHECK(header_size);
    _CHECK(library_version);
    _CHECK(simulator_version);
    _CHECK(amount_of_cells);
    _CHECK(total_amount_of_compartments);
    _CHECK(number_of_steps);
    _CHECK(start_time);
    _CHECK(end_time);
    _CHECK(dt);
    _CHECK(data_unit);
    _CHECK(time_unit);
    _CHECK(mapping_size);
    _CHECK(mapping_name);
    _CHECK(extra_mapping_name);
    _CHECK(extra_mapping_size);
    _CHECK(target_name);

    // else if nothing returned false:
    return true;
#undef _CHECK
}

CellInfo::CellInfo() {
    cell_num = 0;
    master = &cell_num;
    amount_of_compartments = 0;
    data_location = 0;
    extra_mapping_location = 0;
    mapping_location = 0;
}

CellInfo::CellInfo(char* read_buffer, bool is_native, key cmp) {
#define fetch(TYPE, POSITION) *((TYPE*) (read_buffer + POSITION))
    cell_num = fetch(int32_t, SONATA_REPORT_NUMBER_OF_CELL_POSITION);
    amount_of_compartments = fetch(int32_t, SONATA_REPORT_NUMBER_OF_COMPARTMENTS_POSITION);

    data_location = fetch(uint64_t, SONATA_REPORT_DATA_INFO_POSITION);
    extra_mapping_location = fetch(uint64_t, SONATA_REPORT_EXTRA_MAPPING_INFO_POSITION);
    mapping_location = fetch(uint64_t, SONATA_REPORT_MAPPING_INFO_POSITION);

#undef fetch

    if (!is_native) {
        FIX_INT(cell_num);
        FIX_INT(amount_of_compartments);

        DOUBLESWAP(data_location);
        DOUBLESWAP(extra_mapping_location);
        DOUBLESWAP(mapping_location);
    }

    if (cmp == GID)
        master = &cell_num;
    else
        master = &data_location;
}

bool CellInfo::operator==(const CellInfo& c1) const {
#define _CHECK(field)      \
    if (c1.field != field) \
    return false
    _CHECK(cell_num);
    _CHECK(amount_of_compartments);
    _CHECK(data_location);
    _CHECK(extra_mapping_location);
    _CHECK(mapping_location);

    return true;
#undef _CHECK
}

bool CellInfo::operator!=(const CellInfo& c1) const {
    return !(*this == c1);
}

#define COMPARE_(op) \
    (master == &cell_num) ? cell_num op c.cell_num : data_location op c.data_location

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

bool CellInfo::compare_by_offset(const CellInfo& a, const CellInfo& b) {
    return a.data_location < b.data_location;
}

FrameInfo::FrameInfo(char* buffer, char* end_of_buffer, bool is_native, CellInfo::key sort_mode) {
    // create a cellinfo and then insert into itself (rem, this class inherits from set).  The
    // resulting set iterator is then put into
    // a map member for associative lookup gid:location
    for (char* buffer_loc = buffer; buffer_loc < end_of_buffer;
         buffer_loc += SONATA_REPORT_SIZE_CELL_INFO_LENGTH) {
        iterator ins = insert(CellInfo(buffer_loc, is_native, sort_mode)).first;
        gid_lookup.insert(make_pair(ins->get_cell_num(), ins));
    }
}

FrameInfo::FrameInfo(const FrameInfo& rhs)
    : baseTypes::FrameInfo() {
    *this = rhs;
}

FrameInfo& FrameInfo::operator=(const FrameInfo& rhs) {
    insert(rhs.begin(), rhs.end());
    for (iterator cell_item = begin(); cell_item != end(); ++cell_item) {
        gid_lookup.insert(make_pair(cell_item->get_cell_num(), cell_item));
    }

    return *this;
}

FrameInfo& FrameInfo::filter_gids(const vector<CellID>& gids) {
    Lookup new_gid_lookup;

    for (vector<CellID>::const_iterator gid = gids.begin(); gid != gids.end(); gid++) {
        Lookup::iterator it = gid_lookup.find(*gid);
        if (it != gid_lookup.end()) {
            new_gid_lookup.insert(*it);
            gid_lookup.erase(it);
        } else
            throw string("GID doesn't exist in file.");  // TODO: replace by custom std::exception
    }

    // gid_lookup contains the gids that need to be removed and the new_gid_lookup contains the
    // gids that remain.

    for (Lookup::iterator i = gid_lookup.begin(); i != gid_lookup.end(); i++)
        erase(i->second);  // delete the unwanted gids.

    gid_lookup.swap(new_gid_lookup);

    return *this;
}

CompartmentMapping::CompartmentMapping()
    : Base() {}

CompartmentMapping::CompartmentMapping(MappingItem* start, MappingItem* _end)
    : Base(start, _end) {}

CompartmentMapping::CompartmentMapping(MappingItem* buffer, SizeType buf_size)
    : Base(buffer, buffer + buf_size) {}

SectionData::SectionData()
    : Base() {}

SectionData::SectionData(DataItem first_data) {
    push_back(first_data);
}

SectionData::SectionData(DataItem* start, DataItem* _end)
    : Base(start, _end) {}

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

bool SectionData::compare_by_original_location(const Identifier& c1, const Identifier& c2) {
    return c1.get_location() < c2.get_location();
}

CellData::CellData() {
    sect_sort = SID;
}

CellData::CellData(SortingKey sorting) {
    sect_sort = sorting;
}

CellData::CellData(const DataItem* data_buffer, const CompartmentMapping& mapping) {
    sect_sort = SID;
    fill_data(data_buffer, mapping);
}

void CellData::fill_data(const DataItem* data_buffer, const CompartmentMapping& mapping) {
    CompartmentMapping::const_iterator it;
    SectionIndex i;
    for (it = mapping.begin(), i = 0; it != mapping.end(); it++, data_buffer++, i++) {
        // iteration through all compartments inside cell
        DataItem data = *data_buffer;

        iterator seg_match = find(*it);  // find corresponding segment

        if (seg_match == end()) {
            // SID was not found.
            // this means it wasnt created yet.
            // therefore create a new one and insert first data item into it:
            iterator ins =
                insert(make_pair(SectionData::Identifier(*it, i, sect_sort), SectionData(data)))
                    .first;
            // note that it is guaranteed that the above insertion takes place, because this
            // block is only entered
            // if the sid does not exist inside this CellData object yet.
            // if sorting is done by location, this will be unique anyways.
            sid_lookup.insert(make_pair(*it, ins));
        } else
            // SID is already present in CellData. The new data item just needs to be inserted
            seg_match->second.push_back(data);
    }
}

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

ReadManager::ReadManager() {}

ReadManager::ReadManager(const string& file_path)
    : ifstream(file_path.c_str(), ios::in | ios::binary) {
    exceptions(ios::badbit | ios::eofbit | ios::failbit);
    read_header();
}

void ReadManager::open(const string& file_path) {
    exceptions(ios::badbit | ios::eofbit | ios::failbit);

    ifstream::open(file_path.c_str(), ios::in | ios::binary);
    read_header();
}

void ReadManager::read_header() {
    char read_buffer[SONATA_REPORT_HEADER_LENGTH + 1];

    seekg(ios::beg);

    ifstream::read(read_buffer, SONATA_REPORT_HEADER_LENGTH);

    read_buffer[SONATA_REPORT_HEADER_LENGTH] = 0;  // c string terminator
    file_info = read_buffer;

    cell_info_origin = tellg();

    is_native = file_info.source_is_native();

    CellInfo first_cell = retrieve_cell_info(0);

    start_location = first_cell.get_data_location();
    mapping_location = first_cell.get_mapping_location();

    mapping_size_floats = file_info.get_mapping_size() *
                          file_info.get_total_number_of_compartments();
}

CellInfo ReadManager::make_null_info() {
    char fake_buffer[SONATA_REPORT_SIZE_CELL_INFO_LENGTH];

    for (char *i = fake_buffer, *buffer_end = fake_buffer + SONATA_REPORT_SIZE_CELL_INFO_LENGTH;
         i < buffer_end;
         i++)
        *i = 0;

    return CellInfo(fake_buffer);
}

const CellInfo ReadManager::_null_info = ReadManager::make_null_info();

CellInfo ReadManager::retrieve_find_cell(CellID cell_num) {
    int start = (cell_num > file_info.get_number_of_cells())
                    ? file_info.get_number_of_cells() - 1
                    : cell_num - 1;  // in many cases cellNum will be same as cellIndex
    CellInfo seek;

    bool opposite_search = false;

    int i = start;
    do {
        if (i >= (int) file_info.get_number_of_cells() || i < 0) {
            if (opposite_search)
                return make_null_info();
            else {
                opposite_search = true;
                i = start;
            }
        }

        seek = retrieve_cell_info(i);

        if ((int) seek.get_cell_num() < i)
            if (opposite_search)
                i--;
            else
                i++;
        else if (opposite_search)
            i++;
        else
            i--;

    } while (seek.get_cell_num() != cell_num);

    return seek;
}

CellInfo ReadManager::retrieve_cell_info(CellIndex cell_index, CellInfo::key sort_mode) {
    FilePosition ppos = tellg();

    seekg_cell_info(cell_index);

    char read_buffer[SONATA_REPORT_SIZE_CELL_INFO_LENGTH];
    read(read_buffer, SONATA_REPORT_SIZE_CELL_INFO_LENGTH);

    CellInfo req_info(read_buffer, is_native, sort_mode);

    seekg(ppos);

    return req_info;
}

FrameInfo ReadManager::retrieve_all_cell_info(CellInfo::key sort_mode) {
    FilePosition ppos = tellg();

    CellCount amount_of_cells = file_info.get_number_of_cells();
    streamsize bytes_to_read = SONATA_REPORT_SIZE_CELL_INFO_LENGTH * amount_of_cells;
    char* read_buffer = new char[bytes_to_read];
    char* end_of_buffer = read_buffer + bytes_to_read;
    seekg_cell_info(0);

    read(read_buffer, bytes_to_read);

    FrameInfo all_cells(read_buffer, end_of_buffer, is_native, sort_mode);

    delete[] read_buffer;

    seekg(ppos);

    return all_cells;
}

void ReadManager::read_frame(DataItem* data_buffer) {
    read(data_buffer, file_info.get_total_number_of_compartments());
}

void ReadManager::read_frame_mapping(MappingItem* buffer) {
    FilePosition ppos = tellg();

    seekg_mapping_start();

    read(buffer, mapping_size_floats);

    seekg(ppos);
}

FrameParser::FrameOrganiser::FrameOrganiser() {
    sect_sort = SID;
}

FrameParser::FrameOrganiser::FrameOrganiser(SortingKey segment_sorting) {
    sect_sort = segment_sorting;
}

void FrameParser::FrameOrganiser::fill_data(const FrameInfo& info,
                                            const vector<CompartmentMapping>& mapping,
                                            const DataItem* data_buffer) {
    FrameInfo::const_iterator pinf = info.begin();
    vector<CompartmentMapping>::const_iterator pmap = mapping.begin();

    for (const DataItem* ptr = data_buffer; pinf != info.end();
         ptr += pinf->get_amount_of_compartments(), pinf++, pmap++) {
        // iteration through all gids

        iterator ins = insert(make_pair(*pinf, CellData(sect_sort))).first;  // insert CellData
                                                                             // object for current
                                                                             // gid
        ins->second.fill_data(ptr, *pmap);  // feed new CellData object with mapping and data
    }
}

void FrameParser::FrameOrganiser::make_ordered_data_buffer(DataItem* buffer) const {
    for (const_iterator cd = begin(); cd != end(); cd++)  // iterate through cells
        for (CellData::const_iterator sd = cd->second.begin(); sd != cd->second.end();
             sd++)  // iterate through segments
            for (SectionData::const_iterator f = sd->second.begin(); f != sd->second.end();
                 f++, buffer++)  // iterate through compartments
                *buffer = *f;
}

FrameParser::CompartmentIndexing FrameParser::FrameOrganiser::make_ordered_offset_references()
    const {
    CompartmentIndexing gid_v;
    gid_v.reserve(size());

    FrameDataIndex position = 0;

    for (const_iterator gid = begin(); gid != end(); gid++) {  // iterate through gids
        gid_v.push_back(vector<FrameDataIndex>());

        // reserve memory for amount of segments
        gid_v.back().reserve((int) (gid->second.size() * 1.5));

        MappingItem presid = 0;
        for (CellData::const_iterator sid = gid->second.begin(); sid != gid->second.end();
             position += sid->second.size(), sid++) {
            // iterating through available sid's
            while (presid++ != sid->first.get_sid())
                gid_v.back().push_back(0);  // insert dummy variables for unavailable sids (gaps)
            gid_v.back().push_back(position);
        }
    }

    return gid_v;
}

FrameParser::CompartmentCounts FrameParser::FrameOrganiser::make_ordered_seg_sizes() const {
    CompartmentCounts gid_v;
    gid_v.reserve(size());

    for (const_iterator gid = begin(); gid != end(); gid++) {
        gid_v.push_back(vector<SectionDataCount>());
        gid_v.back().reserve((int) (gid->second.size() * 1.5));
        MappingItem presid = 0;
        for (CellData::const_iterator sid = gid->second.begin(); sid != gid->second.end(); sid++) {
            while (presid++ != sid->first.get_sid())
                gid_v.back().push_back(0);  // fill sid gaps with zero sizes
            gid_v.back().push_back(sid->second.size());
        }
    }

    return gid_v;
}

FrameParser::FrameParser() {
    read_buffer = NULL;
    ref_array = NULL;
}

FrameParser::FrameParser(const string& file, bool sort_data)
    : ReadManager(file) {
    if (file_info.get_mapping_size() != 1) {
        cerr << "Target file has incompatible mapping size for FrameParser. FrameParser only works "
                "with mapping size 1."
             << endl;
        exit(1);
    }

    CellInfo::key cinfo_sorting = (sort_data) ? CellInfo::GID : CellInfo::LOCATION;

    FrameInfo cinfo = retrieve_all_cell_info(cinfo_sorting);

    ReadOptimiser::set_default_element_size(sizeof(DataItem));
    offsets = ReadOptimiser(file_info.get_total_number_of_compartments());
    create_references(cinfo, sort_data);

    timestep = 0;  // by default start pointing to timstep 0
    ReadManager::seekg_frame(timestep);
}

FrameParser::FrameParser(const string& file, const vector<CellID>& gid_target, bool sort_data)
    : ReadManager(file) {
    timestep = 0;
    ReadOptimiser::set_default_element_size(sizeof(DataItem));
    retarget(gid_target, sort_data);
}

FrameParser::~FrameParser() {
    if (read_buffer)  // readBuffer may have not been used if sorting was disabled, in which case
                      // it would have been set to NULL
        delete[] read_buffer;
    if (ref_array)  // same applies for refArray
        delete[] ref_array;
}

void FrameParser::retarget(const vector<CellID>& gid_target, bool sort_data) {
    CellInfo::key cinfo_sorting = (sort_data) ? CellInfo::GID : CellInfo::LOCATION;

    FrameInfo cinfo = retrieve_all_cell_info(cinfo_sorting);
    cinfo.filter_gids(gid_target);

    offsets = ReadOptimiser(cinfo, get_frame_size_bytes(), get_start_offset());
    offsets.set_element_size(sizeof(DataItem));

    create_references(cinfo, sort_data);

    seekg(get_start_offset() + offsets.get_first());

    operator=(timestep);
}

void FrameParser::create_references(const FrameInfo& cinfo, bool sort_data) {
    SortingKey s_sorting = (sort_data) ? SID : LOC;

    ref_array = NULL;
    read_buffer = NULL;

    vector<CompartmentMapping> all_mapping;

    DataItem* temp_buffer =
        new DataItem[offsets.get_data_size_elements()];  // allocate memory for readBuffer. This
                                                         // will atleast be needed for creating
                                                         // the FrameOrganiser

    read_frame_mapping(temp_buffer);

    FrameInfo::const_iterator pinf = cinfo.begin();

    for (DataItem *mpptr = temp_buffer, *endptr = mpptr + pinf->get_amount_of_compartments();
         pinf != cinfo.end();
         pinf++, mpptr = endptr, endptr += pinf->get_amount_of_compartments())
        all_mapping.push_back(CompartmentMapping(mpptr, endptr));

    if (sort_data)
        for (FrameDataIndex i = 0; i < offsets.get_data_size_elements(); i++)
            temp_buffer[i] = i;  // makes buffer of indexes { 0 , 1 ,2 , 3 , ... }
    // else it doesnt really matter what FrameOrganiser is fed with.

    org = FrameOrganiser(s_sorting);

    org.fill_data(cinfo, all_mapping, temp_buffer);  // feed the FrameOrganiser with the data.

    if (sort_data) {
        read_buffer = temp_buffer;  // keep allocated memory at tempBuffer

        // sort_data is enabled -> create an array of pointers to readBuffer
        // in which the pointers are sorted by GID and SID of the data the point to.
        // these pointers are to be stored into ref array
        org.make_ordered_data_buffer(read_buffer);  // rearranges indexes to correspond to gid and
                                                    // segment order
        // the memory allocated for readBuffer is now being used for a different purpose:
        // it contains the sorted indexes, which now need to be converted to pointers to
        // readBuffer.

        ref_array = new DataItem*[offsets.get_data_size_elements()];  // allocate memory for
                                                                      // refArray (will contain
                                                                      // pointers for sorting)

        for (float *p_index = read_buffer,
                   **p_ptr = ref_array,
                   *p_end = read_buffer + offsets.get_data_size_elements();
             p_index != p_end;
             p_index++, p_ptr++)
            *p_ptr = &read_buffer[(FrameDataIndex) *p_index];  // converts indexes to pointers and
                                                               // stores them into refArray
    } else
        // data sorting is disabled -> allocated memory at temp Buffer is no longer needed.
        //(data will be read directly into the memory that user has provided)
        delete[] temp_buffer;
}

void FrameParser::buffer_transfer(DataItem* buffer) const {
    DataItem* req_location = buffer;

    if (buffer == read_buffer)  // source and destination is the same - trouble!
        buffer = new DataItem[offsets.get_data_size_elements()];  // need to allocate different
                                                                  // destination

    for (DataItem **source_ptr = ref_array,
                  *target = buffer,
                  *d_end = buffer + offsets.get_data_size_elements();
         target != d_end;
         source_ptr++, target++)
        *target = **source_ptr;

    if (req_location != buffer) {  // buffer has not been transfered where requested!
        // this can now be corrected.
        for (DataItem *origin = buffer,
                      *destination = req_location,
                      *d_end = buffer + offsets.get_data_size_elements();
             origin != d_end;
             origin++, destination++)
            *destination = *origin;  // copy data to requested location
        delete[] buffer;  // buffer was is only a temporary memory allocation - delete it.
    }
}

FrameParser& FrameParser::operator++() {
    timestep++;
    seekg(file_info.get_total_number_of_compartments() * sizeof(DataItem), cur);
    return *this;
}

FrameParser& FrameParser::operator++(int) {
    return this->operator++();
}

FrameParser& FrameParser::operator--() {
    timestep--;
    seekg(-file_info.get_total_number_of_compartments() * sizeof(DataItem), cur);
    return *this;
}

FrameParser& FrameParser::operator--(int) {
    return this->operator--();
}

FrameParser& FrameParser::operator+=(FrameDiff increment) {
    timestep += increment;
    seekg(increment * file_info.get_total_number_of_compartments() * sizeof(DataItem), cur);
    return *this;
}

FrameParser& FrameParser::operator-=(FrameDiff decrement) {
    timestep -= decrement;
    seekg(-decrement * file_info.get_total_number_of_compartments() * sizeof(DataItem), cur);
    return *this;
}

FrameParser& FrameParser::operator=(FrameIndex new_time) {
    timestep = new_time;
    seekg_frame(timestep);
    seekg(offsets.get_first(), cur);
    return *this;
}

Time FrameParser::get_time() const {
    return file_info.get_start_time() + timestep * file_info.get_time_step_size();
}

void FrameParser::read_frame(DataItem* buffer) {
    DataItem* read_ptr;

    if (read_buffer)             // readBuffer will be NULL if sorting is disabled
        read_ptr = read_buffer;  // sorting enabled -> read into readBuffer and transfer data in
                                 // order later
    else
        read_ptr = buffer;  // sorting disabled -> read directly into buffer

    // offsets contains alternating values of data sizes to be read and skipped.
    for (ReadOptimiser::iterator i = offsets.begin(); i != offsets.end(); i++) {
        read(read_ptr, *i);
        read_ptr += *i;
        i++;

        // seems to be confusion as to whether offsets should contain bytelength or element
        // count.  I have it contain element count for now
        // i now points to an offset that indicates amount of data to be skipped
        if (i == offsets.end())  // no more offsets -> i'm done reading
            break;

        if (*i)
            seekg(*i, ios::cur);  // non-zero offset -> perfom the skip
    }

    if (ref_array)                // refArray will be NULL if sorting is disabled
        buffer_transfer(buffer);  // sorting is enabled. Transfer data from readBuffer to buffer
                                  // such that data is in order
}

void FrameParser::read_frame_mapping(MappingItem* buffer) {
    // readFrameMapping should not return without the filepointer being in the same state as
    // before it was called.
    FilePosition before = tellg();  // make a backup for current filepointer

    seekg_mapping_start();  // now move file pointer to mapping location

    size_t prev_element = offsets.get_element_size();  // backup element size.
    offsets.set_element_size(sizeof(MappingItem));     // setup offsets to be reading MappingItems
    // instead of whatever it was reading before.

    seekg(offsets.get_first(), ios::cur);  // and advance to first cell
    read_frame(buffer);  // read the data - this function will also alter the filepointer and
                         // timestep

    offsets.set_element_size(prev_element);  // restore element size.

    seekg(before);  // restore filepointer
}

void FrameParser::read_frame_data(DataItem* buffer) {
    read_frame(buffer);
    timestep++;  // the filepointer has been incremented to the next timestep as data has been
                 // read.
}

size_t FrameParser::ReadOptimiser::default_element_size = 1;

FrameParser::ReadOptimiser::ReadOptimiser() {
    element_size = 0;
    first_offset = 0;
    data_size_elements = 0;
}

FrameParser::ReadOptimiser::ReadOptimiser(std::streamsize frame_size_bytes) {
    element_size = default_element_size;
    first_offset = 0;
    data_size_elements = frame_size_bytes;  // * element_size;
    push_back(data_size_elements);
}

FrameParser::ReadOptimiser::ReadOptimiser(const FrameInfo& cells_to_read,
                                          std::streamsize frame_size_bytes,
                                          FilePosition start_location) {
    element_size = default_element_size;
    // collect data edge points
    multiset<FilePosition> data_edges;

    data_edges.insert(start_location);

    for (FrameInfo::const_iterator i = cells_to_read.begin(); i != cells_to_read.end(); i++) {
        data_edges.insert(i->get_data_location());  // starting point
        data_edges.insert(i->get_data_location() +
                          i->get_amount_of_compartments() * element_size);  // end point
    }

    data_edges.insert(start_location + (FileOffset) frame_size_bytes);

    /*create a vector that contains the offsets.
     * The first offset represents the first bunch of data to be read.
     * The second offset then represents the amound of data to be skipped after the first read.
     * and then the next offset is again the amount of data to be read.
     */

    // compute the offsets differentiating the dataEdges
    list<FileOffset> offsets;
    for (multiset<FilePosition>::iterator current_pos = data_edges.begin(),
                                          prev_pos = current_pos++;
         current_pos != data_edges.end();
         prev_pos = current_pos, current_pos++)
        offsets.push_back(*current_pos - *prev_pos);

    // the last offsets should move the file pointer beyond the frame to the start location on
    // the next frame
    // therefore the first and last offsets are special -> take them out

    first_offset = offsets.front();
    FileOffset last_offset = offsets.back();

    offsets.pop_front();  // remove first offset from list
    offsets.pop_back();   // remove last offset from list

    // now the first an last offset are read dimensions. This guarantees that a skip has always
    // a precessor and successor
    // this is convenient in the next step, as a skip is always on an even location aswell as
    // offsets.end()

    list<FileOffset>::iterator it;
    bool toggle_on_read;

    data_size_elements = 0;
    for (toggle_on_read = true, it = offsets.begin(); it != offsets.end();
         toggle_on_read = !toggle_on_read, it++) {
        if (toggle_on_read) {
            data_size_elements += (*it /= element_size);
            // handle read
            // read specifications should be in units of floats, but skips should be in units of
            // bytes.
            // also size of data to be read per frame needs to established
        } else
            // handle skip - There should be no zero offset inbetween
            if (*it == 0) {
            // remove this zero skip by combining prev. and next read.
            list<FileOffset>::iterator prev_read = it;
            list<FileOffset>::iterator next_read = it;
            prev_read--;
            next_read++;

            data_size_elements += (*next_read /= element_size);  // next read was still in units of
                                                                 // bytes
            *prev_read += *next_read;                            // combine nextRead with prevRead
            offsets.erase(it);
            offsets.erase(next_read);
            it = prev_read;         // point to prevRead now
            toggle_on_read = true;  // and mark that toggle flag that current iteration deal with
                                    // read.
        }
    }

    // now optimise access performance by converting the offsets list to a vector.

    vector<FileOffset> voffsets(offsets.begin(), offsets.end());
    swap(voffsets);

    // compute the finalSkip, which brings the filePointer back to the starting location, but of
    // the nextFrame
    FileOffset final_skip = last_offset + first_offset;

    if (final_skip)
        push_back(final_skip);
}

void FrameParser::ReadOptimiser::update_offsets(size_t new_size) {
    data_size_elements = new_size * data_size_elements / element_size;
    for (iterator i = begin(); i != end(); i++) {  // iterate through all read offsets
        // first reduce the read offsets to byte numbers, then convert it to the new element
        // size
        *i = new_size * (*i) / element_size;
        i++;
        if (i == end())
            break;
    }
}

}  // namespace binary_reader
}  // namespace bbp
