#pragma once

// Included after counter and the transient barcode_entry are defined.
// Resident entries have a fixed word width; temporary queries may widen.
template <class Word> struct basic_barcode_entry {
    static_assert(std::is_same_v<Word, uint64_t> || std::is_same_v<Word, uint128_t>,
                  "Resident barcodes require an unsigned 64-bit or 128-bit word");
    Word bits{};
    mutable counter count{};
    uint16_t length{};
    basic_barcode_entry() = default;
    static Word word_bits(const int64_seq& sequence) {
        if constexpr (std::is_same_v<Word, uint64_t>)
            return sequence.as_uint64();
        else
            return sequence.as_uint128();
    }
    explicit basic_barcode_entry(const barcode_entry& entry)
        : bits(word_bits(entry.barcode)), count(entry.count), length(entry.barcode.length) {}
    explicit basic_barcode_entry(const packed_barcode<Word>& key)
        : bits(key.bits), length(key.length) {}
    packed_barcode<Word> key() const noexcept { return {bits, length}; }
    int64_seq sequence() const { return int64_seq::from_word(bits, length); }
    bool operator==(const basic_barcode_entry& other) const noexcept {
        return bits == other.bits && length == other.length;
    }
    bool is_valid() const noexcept { return length > 0 && length <= sizeof(Word) * 4; }
};
template <class Word> struct basic_barcode_entry_hash {
    size_t operator()(const basic_barcode_entry<Word>& entry) const noexcept {
        return std::hash<packed_barcode<Word>>{}(entry.key());
    }
};

class barcode_entry_handle {
    const void* pointer_ = nullptr;
    bool wide_ = false;

  public:
    barcode_entry_handle() = default;
    template <class Word>
    explicit barcode_entry_handle(const basic_barcode_entry<Word>* entry)
        : pointer_(entry), wide_(std::is_same_v<Word, uint128_t>) {}
    explicit operator bool() const noexcept { return pointer_ != nullptr; }
    template <class F> decltype(auto) visit(F&& fn) const {
        if (wide_)
            return fn(*static_cast<const basic_barcode_entry<uint128_t>*>(pointer_));
        return fn(*static_cast<const basic_barcode_entry<uint64_t>*>(pointer_));
    }
    int64_seq barcode() const {
        if (!pointer_)
            return {};
        return visit([](const auto& entry) { return entry.sequence(); });
    }
    counter* counts() const {
        if (!pointer_)
            return nullptr;
        return visit([](const auto& entry) { return &entry.count; });
    }
    barcode_entry snapshot() const {
        barcode_entry result;
        result.barcode = barcode();
        if (const auto* c = counts())
            result.count = *c;
        return result;
    }
    bool operator==(const barcode_entry_handle& other) const noexcept {
        return pointer_ == other.pointer_ && wide_ == other.wide_;
    }
    bool operator!=(const barcode_entry_handle& other) const noexcept { return !(*this == other); }
    size_t hash() const noexcept { return std::hash<const void*>{}(pointer_) ^ static_cast<size_t>(wide_); }
};
struct barcode_entry_handle_hash {
    size_t operator()(const barcode_entry_handle& entry) const noexcept { return entry.hash(); }
};

namespace std {
template <> struct hash<barcode_entry_handle> {
    size_t operator()(const barcode_entry_handle& entry) const noexcept { return entry.hash(); }
};
} // namespace std

namespace barcode_storage_detail {
template <class Narrow, class Wide> class table_pair {
    Narrow narrow_;
    std::unique_ptr<Wide> wide_owner_;
    std::atomic<Wide*> wide_{nullptr};
    mutable std::mutex initialization_mutex_;

  public:
    table_pair() = default;
    table_pair(const table_pair& other) : narrow_(other.narrow_) {
        if (const auto* source = other.wide()) {
            wide_owner_ = std::make_unique<Wide>(*source);
            wide_.store(wide_owner_.get(), std::memory_order_release);
        }
    }
    table_pair(table_pair&& other)
        : narrow_(std::move(other.narrow_)), wide_owner_(std::move(other.wide_owner_)) {
        wide_.store(wide_owner_.get(), std::memory_order_release);
        other.wide_.store(nullptr, std::memory_order_release);
    }
    table_pair& operator=(const table_pair& other) {
        if (this != &other) {
            table_pair copy(other);
            *this = std::move(copy);
        }
        return *this;
    }
    table_pair& operator=(table_pair&& other) {
        if (this != &other) {
            narrow_ = std::move(other.narrow_);
            wide_owner_ = std::move(other.wide_owner_);
            wide_.store(wide_owner_.get(), std::memory_order_release);
            other.wide_.store(nullptr, std::memory_order_release);
        }
        return *this;
    }
    Narrow& narrow() noexcept { return narrow_; }
    const Narrow& narrow() const noexcept { return narrow_; }
    Wide* wide() noexcept { return wide_.load(std::memory_order_acquire); }
    const Wide* wide() const noexcept { return wide_.load(std::memory_order_acquire); }
    Wide& ensure_wide() {
        if (auto* existing = wide())
            return *existing;
        std::lock_guard<std::mutex> lock(initialization_mutex_);
        if (!wide_owner_) {
            wide_owner_ = std::make_unique<Wide>();
            wide_.store(wide_owner_.get(), std::memory_order_release);
        }
        return *wide_owner_;
    }
    size_t size() const {
        const auto* w = wide();
        return narrow_.size() + (w ? w->size() : 0);
    }
    void clear() {
        // Clearing, moving and destruction require quiescent readers.
        narrow_.clear();
        wide_.store(nullptr, std::memory_order_release);
        wide_owner_.reset();
    }
    template <class F> void for_each(F&& fn) const {
        narrow_.for_each(fn);
        if (const auto* w = wide())
            w->for_each(fn);
    }
};
template <class Word> packed_barcode<Word> pack_key(const int64_seq& sequence) {
    if constexpr (std::is_same_v<Word, uint64_t>)
        return {sequence.as_uint64(), sequence.length};
    else
        return {sequence.as_uint128(), sequence.length};
}
inline const int64_seq& key_of(const int64_seq& sequence) {
    return sequence;
}
inline const int64_seq& key_of(const barcode_entry& entry) {
    return entry.barcode;
}
inline int64_seq key_of(const barcode_entry_handle& entry) {
    return entry.barcode();
}
template <class Word> int64_seq key_of(const packed_barcode<Word>& sequence) {
    return int64_seq(sequence);
}
template <class Word> int64_seq key_of(const basic_barcode_entry<Word>& entry) {
    return entry.sequence();
}
template <class T> auto key_of(const T* pointer) -> decltype(key_of(*pointer)) {
    return key_of(*pointer);
}
template <class T, class = void> struct is_key_input : std::false_type {};
template <class T>
struct is_key_input<T, std::void_t<decltype(key_of(std::declval<const T&>()))>> : std::true_type {};
template <class T, class F> void for_each_key(const T& input, F&& fn) {
    if constexpr (is_key_input<T>::value)
        fn(key_of(input));
    else
        for (const auto& element : input)
            for_each_key(element, fn);
}
using summary_row = std::tuple<std::string, int, int, int>;
inline summary_row summarize(const int64_seq& sequence, const counter* counts) {
    const counter zero;
    const auto& c = counts ? *counts : zero;
    return {sequence.bits_to_sequence(), c.load(raw), c.load(total), c.load(corrected)};
}
inline void write_header(std::ostream& out) {
    out << "class_id,true_barcode,raw_count,total_count,corrected_count\n";
}
inline void write_row(std::ostream& out, const std::string& class_id, const summary_row& row) {
    out << class_id << ',' << std::get<0>(row) << ',' << std::get<1>(row) << ',' << std::get<2>(row) << ','
        << std::get<3>(row) << '\n';
}
template <class Word>
using flat_entries = phmap::parallel_flat_hash_set<basic_barcode_entry<Word>, basic_barcode_entry_hash<Word>,
                                                   std::equal_to<basic_barcode_entry<Word>>>;
template <class Word>
using node_entries = phmap::parallel_node_hash_set<basic_barcode_entry<Word>, basic_barcode_entry_hash<Word>,
                                                   std::equal_to<basic_barcode_entry<Word>>,
                                                   std::allocator<basic_barcode_entry<Word>>, 6, std::mutex>;
using handle_set = phmap::flat_hash_set<barcode_entry_handle, barcode_entry_handle_hash>;
template <class Word>
using alias_map = phmap::parallel_node_hash_map<
    packed_barcode<Word>, handle_set, std::hash<packed_barcode<Word>>, std::equal_to<packed_barcode<Word>>,
    std::allocator<std::pair<const packed_barcode<Word>, handle_set>>, 6, std::mutex>;
class handle_iterator {
    std::shared_ptr<std::vector<barcode_entry_handle>> entries_;
    size_t position_ = 0;

  public:
    handle_iterator() = default;
    handle_iterator(std::shared_ptr<std::vector<barcode_entry_handle>> entries, size_t position)
        : entries_(std::move(entries)), position_(position) {}
    const barcode_entry_handle& operator*() const { return (*entries_)[position_]; }
    const barcode_entry_handle* operator->() const { return &(*entries_)[position_]; }
    handle_iterator& operator++() {
        ++position_;
        return *this;
    }
    bool operator==(const handle_iterator& other) const {
        return entries_ == other.entries_ && position_ == other.position_;
    }
    bool operator!=(const handle_iterator& other) const { return !(*this == other); }
};
inline auto handle_range(std::vector<barcode_entry_handle> entries) {
    auto storage = std::make_shared<std::vector<barcode_entry_handle>>(std::move(entries));
    return std::make_pair(handle_iterator(storage, 0), handle_iterator(storage, storage->size()));
}
} // namespace barcode_storage_detail

struct bc_flat_set {
    using key = int64_seq;
    using value = barcode_entry;

  private:
    // Structural mutation is quiescent, as for the resident NullMutex tables.
    // Erasure may retain old bits: this mask is a conservative length filter.
    uint64_t length_mask_ = 0;
    barcode_storage_detail::table_pair<barcode_storage_detail::flat_entries<uint64_t>,
                                       barcode_storage_detail::flat_entries<uint128_t>>
        entries_;
    barcode_entry_handle find_handle(const key& sequence) const {
        if (!sequence.is_valid())
            return {};
        barcode_entry_handle result;
        if (sequence.is_wide()) {
            if (const auto* table = entries_.wide()) {
                basic_barcode_entry<uint128_t> probe(
                    barcode_storage_detail::pack_key<uint128_t>(sequence));
                table->if_contains(probe, [&](const auto& entry) { result = barcode_entry_handle(&entry); });
            }
        } else {
            basic_barcode_entry<uint64_t> probe(
                barcode_storage_detail::pack_key<uint64_t>(sequence));
            entries_.narrow().if_contains(probe,
                                          [&](const auto& entry) { result = barcode_entry_handle(&entry); });
        }
        return result;
    }

  public:
    size_t size() const { return entries_.size(); }
    size_t association_size() const { return size(); }
    size_t unique_val_size() const { return size(); }
    size_t narrow_size() const { return entries_.narrow().size(); }
    size_t wide_size() const {
        const auto* w = entries_.wide();
        return w ? w->size() : 0;
    }
    bool has_wide_storage() const { return entries_.wide() != nullptr; }
    bool empty() const { return size() == 0; }
    uint64_t length_mask() const noexcept { return length_mask_; }
    void clear() { entries_.clear(); length_mask_ = 0; }
    template <class T> bool check_wl_for(const T& input) const {
        bool found = false;
        barcode_storage_detail::for_each_key(input, [&](const key& sequence) {
            if (!found)
                found = bool(find_handle(sequence));
        });
        return found;
    }
    bool has_identity_mapping(const key& sequence) const { return bool(find_handle(sequence)); }
    template <class T> std::unordered_set<key> return_matching_barcodes(const T& input) const {
        std::unordered_set<key> result;
        barcode_storage_detail::for_each_key(input, [&](const key& sequence) {
            if (find_handle(sequence))
                result.insert(sequence);
        });
        return result;
    }
    template <class T> std::unordered_set<key> return_putative_correct_bcs(const T& input) const {
        return return_matching_barcodes(input);
    }
    template <class T> void insert_bc_entry(const T&, const value& correct) {
        if (!correct.is_valid())
            return;
        // Record the stored identity, not the possibly different observed key.
        length_mask_ |= uint64_t{1} << (correct.barcode.length - 1);
        if (correct.barcode.is_wide())
            entries_.ensure_wide().emplace(basic_barcode_entry<uint128_t>(correct));
        else
            entries_.narrow().emplace(basic_barcode_entry<uint64_t>(correct));
    }
    template <class T> void insert_bc_entry(const T& input, value&& correct) {
        insert_bc_entry(input, static_cast<const value&>(correct));
    }
    template <class T> void insert_bc_entry(const T& input) {
        value correct;
        correct.barcode = barcode_storage_detail::key_of(input);
        insert_bc_entry(input, correct);
    }
    template <class K, class V> void emplace(K&& observed, V&& correct) {
        insert_bc_entry(observed, std::forward<V>(correct));
    }
    template <class T> size_t erase(const T& input) {
        const auto& sequence = barcode_storage_detail::key_of(input);
        if (!sequence.is_valid())
            return 0;
        if (sequence.is_wide()) {
            auto* table = entries_.wide();
            if (!table)
                return 0;
            basic_barcode_entry<uint128_t> probe(
                barcode_storage_detail::pack_key<uint128_t>(sequence));
            return table->erase(probe);
        }
        basic_barcode_entry<uint64_t> probe(
            barcode_storage_detail::pack_key<uint64_t>(sequence));
        return entries_.narrow().erase(probe);
    }
    template <class T> void remove_bc_entry(const T& input) { erase(input); }
    template <class T> counter* find_counter(const T& input) const {
        return find_handle(barcode_storage_detail::key_of(input)).counts();
    }
    template <class T> counter get_all_bc_counts(const T& input) const {
        const auto* counts = find_counter(input);
        return counts ? *counts : counter{};
    }
    template <class T> void update_bc_count(const T& input, barcode_counts slot = total) const {
        if (auto* counts = find_counter(input))
            counts->increment(slot);
    }
    template <class T> void subtract_bc_count(const T& input, barcode_counts slot = total) const {
        if (auto* counts = find_counter(input))
            counts->subtract(slot);
    }
    template <class T> int get_bc_count(const T& input, barcode_counts slot = total) const {
        const auto* counts = find_counter(input);
        return counts ? counts->load(slot) : 0;
    }
    template <class T> void set_bc_count(const T& input, int count, barcode_counts slot = total) const {
        if (auto* counts = find_counter(input))
            counts->counts[slot].store(count, std::memory_order_relaxed);
    }
    template <class F> void for_each_entry(F&& fn) const { entries_.for_each(fn); }
    std::vector<barcode_entry_handle> get_unique_entries() const {
        std::vector<barcode_entry_handle> result;
        result.reserve(size());
        for_each_entry([&](const auto& entry) { result.emplace_back(&entry); });
        return result;
    }
    auto equal_range(const key& sequence) const {
        std::vector<barcode_entry_handle> result;
        if (const auto entry = find_handle(sequence))
            result.push_back(entry);
        return barcode_storage_detail::handle_range(std::move(result));
    }
    std::vector<barcode_storage_detail::summary_row>
    summarize_counts(const std::unordered_set<key>* filter_keys = nullptr) const {
        std::vector<barcode_storage_detail::summary_row> rows;
        rows.reserve(filter_keys ? filter_keys->size() : size());
        if (filter_keys) {
            for (const auto& sequence : *filter_keys)
                if (const auto entry = find_handle(sequence))
                    rows.push_back(barcode_storage_detail::summarize(sequence, entry.counts()));
        } else
            for_each_entry([&](const auto& entry) {
                rows.push_back(barcode_storage_detail::summarize(entry.sequence(), &entry.count));
            });
        return rows;
    }
    void write_wl_summary(std::ostream& out, const std::string class_id,
                          const std::unordered_set<key>* filter_keys = nullptr,
                          bool write_header = true) const {
        if (write_header)
            barcode_storage_detail::write_header(out);
        if (filter_keys) {
            for (const auto& sequence : *filter_keys)
                if (const auto entry = find_handle(sequence))
                    barcode_storage_detail::write_row(
                        out, class_id, barcode_storage_detail::summarize(sequence, entry.counts()));
        } else
            for_each_entry([&](const auto& entry) {
                barcode_storage_detail::write_row(
                    out, class_id, barcode_storage_detail::summarize(entry.sequence(), &entry.count));
            });
    }
    size_t unique_storage_bytes() const {
        size_t bytes =
            sizeof(*this) + entries_.narrow().capacity() * (sizeof(basic_barcode_entry<uint64_t>) + 1);
        if (const auto* w = entries_.wide())
            bytes += sizeof(*w) + w->capacity() * (sizeof(basic_barcode_entry<uint128_t>) + 1);
        return bytes;
    }
    size_t association_storage_bytes() const { return 0; }
    size_t validate_association_keys() const noexcept { return 0; }
    size_t validate_association_values() const noexcept { return 0; }
};

template <class Key, class Value> struct bc_multimap {
    static_assert(std::is_same_v<Key, int64_seq> && std::is_same_v<Value, barcode_entry>,
                  "Barcode maps use typed resident entries");
    using key = Key;
    using value = Value;
    using handle_set = barcode_storage_detail::handle_set;

  private:
    barcode_storage_detail::table_pair<barcode_storage_detail::node_entries<uint64_t>,
                                       barcode_storage_detail::node_entries<uint128_t>>
        canonical_;
    barcode_storage_detail::table_pair<barcode_storage_detail::alias_map<uint64_t>,
                                       barcode_storage_detail::alias_map<uint128_t>>
        associations_;
    template <class F> bool with_associations(const key& sequence, F&& fn) const {
        if (!sequence.is_valid())
            return false;
        if (sequence.is_wide()) {
            const auto* table = associations_.wide();
            return table && table->if_contains(barcode_storage_detail::pack_key<uint128_t>(sequence),
                                               [&](const auto& pair) { fn(pair.second); });
        }
        return associations_.narrow().if_contains(barcode_storage_detail::pack_key<uint64_t>(sequence),
                                                  [&](const auto& pair) { fn(pair.second); });
    }
    barcode_entry_handle insert_canonical(const value& correct) {
        barcode_entry_handle result;
        auto insert = [&](auto& table, const auto& candidate) {
            table.lazy_emplace_l(
                candidate, [&](const auto& entry) { result = barcode_entry_handle(&entry); },
                [&](const auto& ctor) {
                    auto* slot = ctor.slot();
                    ctor(candidate);
                    result = barcode_entry_handle(*slot);
                });
        };
        if (correct.barcode.is_wide())
            insert(canonical_.ensure_wide(), basic_barcode_entry<uint128_t>(correct));
        else
            insert(canonical_.narrow(), basic_barcode_entry<uint64_t>(correct));
        return result;
    }
    void add_association(const key& observed, barcode_entry_handle correct) {
        auto insert = [&](auto& table, const auto& packed) {
            table.lazy_emplace_l(
                packed, [&](auto& pair) { pair.second.insert(correct); },
                [&](const auto& ctor) {
                    handle_set values;
                    values.insert(correct);
                    ctor(packed, std::move(values));
                });
        };
        if (observed.is_wide())
            insert(associations_.ensure_wide(), barcode_storage_detail::pack_key<uint128_t>(observed));
        else
            insert(associations_.narrow(), barcode_storage_detail::pack_key<uint64_t>(observed));
    }
    barcode_entry_handle first_handle(const key& sequence) const {
        barcode_entry_handle result;
        with_associations(sequence, [&](const auto& values) {
            if (!values.empty())
                result = *values.begin();
        });
        return result;
    }
    void copy_from(const bc_multimap& other) {
        other.for_each_entry([&](const auto& entry) {
            value snapshot;
            snapshot.barcode = entry.sequence();
            snapshot.count = entry.count;
            insert_canonical(snapshot);
        });
        other.for_each_association([&](const auto& observed, const auto& values) {
            const key sequence(observed);
            for (const auto handle : values)
                add_association(sequence, insert_canonical(handle.snapshot()));
        });
    }

  public:
    bc_multimap() = default;
    bc_multimap(const bc_multimap& other) { copy_from(other); }
    bc_multimap(bc_multimap&&) = default;
    bc_multimap& operator=(const bc_multimap& other) {
        if (this != &other) {
            bc_multimap copy(other);
            *this = std::move(copy);
        }
        return *this;
    }
    bc_multimap& operator=(bc_multimap&&) = default;
    size_t size() const { return associations_.size(); }
    size_t association_size() const { return size(); }
    size_t unique_val_size() const { return canonical_.size(); }
    size_t narrow_size() const { return canonical_.narrow().size(); }
    size_t wide_size() const {
        const auto* w = canonical_.wide();
        return w ? w->size() : 0;
    }
    size_t narrow_association_size() const { return associations_.narrow().size(); }
    size_t wide_association_size() const {
        const auto* w = associations_.wide();
        return w ? w->size() : 0;
    }
    bool has_wide_storage() const { return canonical_.wide() || associations_.wide(); }
    bool empty() const { return size() == 0; }
    void clear_associations() { associations_.clear(); }
    void clear() {
        associations_.clear();
        canonical_.clear();
    }
    template <class T> void insert_bc_entry(const T& observed, const value& correct) {
        const key observed_sequence = barcode_storage_detail::key_of(observed);
        if (!observed_sequence.is_valid() || !correct.is_valid())
            return;
        add_association(observed_sequence, insert_canonical(correct));
    }
    template <class T> void insert_bc_entry(const T& observed, value&& correct) {
        insert_bc_entry(observed, static_cast<const value&>(correct));
    }
    template <class T> void insert_bc_entry(const T& observed) {
        value correct;
        correct.barcode = barcode_storage_detail::key_of(observed);
        insert_bc_entry(observed, correct);
    }
    template <class T, class... Args> void emplace_bc_entry(const T& observed, Args&&... args) {
        value correct{std::forward<Args>(args)...};
        insert_bc_entry(observed, correct);
    }
    template <class K, class V> void emplace(K&& observed, V&& correct) {
        insert_bc_entry(observed, std::forward<V>(correct));
    }
    template <class K, class V> void insert(const std::pair<K, V>& pair) {
        insert_bc_entry(pair.first, pair.second);
    }
    template <class T> size_t erase(const T& input) {
        const auto& sequence = barcode_storage_detail::key_of(input);
        if (!sequence.is_valid())
            return 0;
        if (sequence.is_wide()) {
            auto* table = associations_.wide();
            return table ? table->erase(barcode_storage_detail::pack_key<uint128_t>(sequence)) : 0;
        }
        return associations_.narrow().erase(barcode_storage_detail::pack_key<uint64_t>(sequence));
    }
    template <class T> void remove_bc_entry(const T& input) { erase(input); }
    template <class T> bool check_wl_for(const T& input) const {
        bool found = false;
        barcode_storage_detail::for_each_key(input, [&](const key& sequence) {
            if (!found)
                found = with_associations(sequence, [](const auto&) {});
        });
        return found;
    }
    bool has_identity_mapping(const key& sequence) const {
        bool identity = false;
        with_associations(sequence, [&](const auto& values) {
            for (const auto handle : values)
                if (handle && handle.barcode() == sequence) {
                    identity = true;
                    break;
                }
        });
        return identity;
    }
    size_t count(const key& sequence) const {
        size_t result = 0;
        with_associations(sequence, [&](const auto& values) { result = values.size(); });
        return result;
    }
    auto equal_range(const key& sequence) const {
        std::vector<barcode_entry_handle> result;
        with_associations(sequence, [&](const auto& values) { result.assign(values.begin(), values.end()); });
        return barcode_storage_detail::handle_range(std::move(result));
    }
    template <class T> std::unordered_set<key> return_matching_barcodes(const T& input) const {
        std::unordered_set<key> result;
        barcode_storage_detail::for_each_key(input, [&](const key& sequence) {
            if (check_wl_for(sequence))
                result.insert(sequence);
        });
        return result;
    }
    template <class T> std::unordered_set<key> return_putative_correct_bcs(const T& input) const {
        std::unordered_set<key> result;
        barcode_storage_detail::for_each_key(input, [&](const key& sequence) {
            with_associations(sequence, [&](const auto& values) {
                for (const auto handle : values)
                    if (handle)
                        result.insert(handle.barcode());
            });
        });
        return result;
    }
    template <class T> counter* find_counter(const T& input) const {
        return first_handle(barcode_storage_detail::key_of(input)).counts();
    }
    template <class T> counter get_all_bc_counts(const T& input) const {
        const auto* counts = find_counter(input);
        return counts ? *counts : counter{};
    }
    template <class T> void update_bc_count(const T& input, barcode_counts slot = total) const {
        with_associations(barcode_storage_detail::key_of(input), [&](const auto& values) {
            for (const auto handle : values)
                if (auto* counts = handle.counts())
                    counts->increment(slot);
        });
    }
    template <class T> void subtract_bc_count(const T& input, barcode_counts slot = total) const {
        with_associations(barcode_storage_detail::key_of(input), [&](const auto& values) {
            for (const auto handle : values)
                if (auto* counts = handle.counts())
                    counts->subtract(slot);
        });
    }
    template <class T> int get_bc_count(const T& input, barcode_counts slot = total) const {
        const auto* counts = find_counter(input);
        return counts ? counts->load(slot) : 0;
    }
    template <class T> void set_bc_count(const T& input, int count, barcode_counts slot = total) const {
        if (auto* counts = find_counter(input))
            counts->counts[slot].store(count, std::memory_order_relaxed);
    }
    template <class F> void for_each_entry(F&& fn) const { canonical_.for_each(fn); }
    template <class F> void for_each_association(F&& fn) const {
        associations_.for_each([&](const auto& pair) { fn(pair.first, pair.second); });
    }
    std::vector<barcode_entry_handle> get_unique_entries() const {
        std::vector<barcode_entry_handle> result;
        result.reserve(unique_val_size());
        for_each_entry([&](const auto& entry) { result.emplace_back(&entry); });
        return result;
    }
    std::vector<barcode_storage_detail::summary_row>
    summarize_counts(const std::unordered_set<key>* filter_keys = nullptr) const {
        std::vector<barcode_storage_detail::summary_row> rows;
        rows.reserve(filter_keys ? filter_keys->size() : size());
        if (filter_keys) {
            for (const auto& sequence : *filter_keys)
                rows.push_back(barcode_storage_detail::summarize(sequence, find_counter(sequence)));
        } else
            for_each_association([&](const auto& observed, const auto& values) {
                const auto* counts = values.empty() ? nullptr : values.begin()->counts();
                rows.push_back(barcode_storage_detail::summarize(key(observed), counts));
            });
        return rows;
    }
    void write_wl_summary(std::ostream& out, const std::string class_id,
                          const std::unordered_set<key>* filter_keys = nullptr,
                          bool write_header = true) const {
        if (write_header)
            barcode_storage_detail::write_header(out);
        if (filter_keys) {
            for (const auto& sequence : *filter_keys)
                barcode_storage_detail::write_row(
                    out, class_id, barcode_storage_detail::summarize(sequence, find_counter(sequence)));
        } else
            for_each_association([&](const auto& observed, const auto& values) {
                const auto* counts = values.empty() ? nullptr : values.begin()->counts();
                barcode_storage_detail::write_row(out, class_id,
                                                  barcode_storage_detail::summarize(key(observed), counts));
            });
    }
    void write_wl_summary(const std::string& path = "",
                          const std::unordered_set<key>* filter_keys = nullptr) const {
        std::ofstream out(path);
        if (!out)
            throw std::runtime_error("Failed to open output file: " + path);
        write_wl_summary(out, "", filter_keys);
    }
    size_t unique_storage_bytes() const {
        const auto estimate = [](const auto& table) {
            using entry_type = typename std::decay_t<decltype(table)>::value_type;
            return sizeof(table) + table.capacity() * (sizeof(void*) + 1) +
                   table.size() * (sizeof(entry_type) + 2 * sizeof(void*));
        };
        size_t bytes = estimate(canonical_.narrow());
        if (const auto* w = canonical_.wide())
            bytes += estimate(*w);
        return bytes;
    }
    size_t association_storage_bytes() const {
        const auto estimate = [](const auto& table) {
            using entry_type = typename std::decay_t<decltype(table)>::value_type;
            size_t bytes = sizeof(table) + table.capacity() * (sizeof(void*) + 1) +
                           table.size() * (sizeof(entry_type) + 2 * sizeof(void*));
            table.for_each([&](const auto& pair) {
                bytes += pair.second.capacity() * (sizeof(barcode_entry_handle) + 1);
            });
            return bytes;
        };
        size_t bytes = estimate(associations_.narrow());
        if (const auto* w = associations_.wide())
            bytes += estimate(*w);
        return bytes;
    }
    size_t validate_association_keys() const noexcept {
        size_t invalid = 0;
        try {
            for_each_association([&](const auto& observed, const auto&) {
                if (!observed.length || observed.length > sizeof(observed.bits) * 4)
                    ++invalid;
            });
        } catch (...) {
            return SIZE_MAX;
        }
        return invalid;
    }
    size_t validate_association_values() const noexcept {
        size_t nulls = 0;
        try {
            for_each_association([&](const auto&, const auto& values) {
                for (const auto handle : values)
                    if (!handle)
                        ++nulls;
            });
        } catch (...) {
            return SIZE_MAX;
        }
        return nulls;
    }
};
