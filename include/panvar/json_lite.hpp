#pragma once

// A small JSON reader and string escaper, enough for panvar's own manifests: objects, arrays, strings
// (with \" \\ \/ \b \f \n \r \t and \u escapes below 0x80), integers, true/false/null. Numbers are kept as
// their text so 64-bit integers stay exact.

#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar::json {

struct Value {
    enum class Kind { Null, Bool, Number, String, Array, Object } kind = Kind::Null;
    bool b = false;
    std::string text;                              // a string's value, or a number's digits
    std::vector<Value> items;                      // an array
    std::map<std::string, Value> fields;           // an object

    bool has(const std::string& key) const { return kind == Kind::Object && fields.count(key) != 0; }
    const Value& at(const std::string& key) const {
        const auto it = fields.find(key);
        if (kind != Kind::Object || it == fields.end()) throw std::runtime_error("json: missing key " + key);
        return it->second;
    }
    std::uint64_t u64() const {
        if (kind != Kind::Number || text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
            throw std::runtime_error("json: expected an unsigned integer");
        return std::stoull(text);
    }
    const std::string& str() const {
        if (kind != Kind::String) throw std::runtime_error("json: expected a string");
        return text;
    }
};

class Parser {
public:
    explicit Parser(const std::string& s) : s_(s) {}
    Value parse() {
        Value v = value();
        ws();
        if (i_ != s_.size()) fail("trailing characters");
        return v;
    }

private:
    [[noreturn]] void fail(const std::string& what) { throw std::runtime_error("json: " + what + " at offset " + std::to_string(i_)); }
    void ws() { while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\n' || s_[i_] == '\r' || s_[i_] == '\t')) ++i_; }
    bool eat(char c) { ws(); if (i_ < s_.size() && s_[i_] == c) { ++i_; return true; } return false; }
    void expect(char c) { if (!eat(c)) fail(std::string("expected '") + c + "'"); }
    bool word(const char* w) {
        const std::string t(w);
        if (s_.compare(i_, t.size(), t) == 0) { i_ += t.size(); return true; }
        return false;
    }
    std::string string() {
        expect('"');
        std::string out;
        while (true) {
            if (i_ >= s_.size()) fail("unterminated string");
            const char c = s_[i_++];
            if (c == '"') return out;
            if (c != '\\') { out += c; continue; }
            if (i_ >= s_.size()) fail("unterminated escape");
            const char e = s_[i_++];
            switch (e) {
                case '"': out += '"'; break;
                case '\\': out += '\\'; break;
                case '/': out += '/'; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'n': out += '\n'; break;
                case 'r': out += '\r'; break;
                case 't': out += '\t'; break;
                case 'u': {
                    if (i_ + 4 > s_.size()) fail("short \\u escape");
                    const unsigned long cp = std::stoul(s_.substr(i_, 4), nullptr, 16);
                    if (cp >= 0x80) fail("\\u escape above 0x7f");
                    out += static_cast<char>(cp);
                    i_ += 4;
                    break;
                }
                default: fail("bad escape");
            }
        }
    }
    Value value() {
        ws();
        if (i_ >= s_.size()) fail("unexpected end");
        Value v;
        const char c = s_[i_];
        if (c == '{') {
            ++i_;
            v.kind = Value::Kind::Object;
            if (eat('}')) return v;
            do {
                ws();
                std::string key = string();
                expect(':');
                v.fields[key] = value();
            } while (eat(','));
            expect('}');
        } else if (c == '[') {
            ++i_;
            v.kind = Value::Kind::Array;
            if (eat(']')) return v;
            do v.items.push_back(value()); while (eat(','));
            expect(']');
        } else if (c == '"') {
            v.kind = Value::Kind::String;
            v.text = string();
        } else if (word("true")) { v.kind = Value::Kind::Bool; v.b = true; }
        else if (word("false")) { v.kind = Value::Kind::Bool; v.b = false; }
        else if (word("null")) { v.kind = Value::Kind::Null; }
        else if (c == '-' || (c >= '0' && c <= '9')) {
            v.kind = Value::Kind::Number;
            const std::size_t b = i_;
            ++i_;
            while (i_ < s_.size() && ((s_[i_] >= '0' && s_[i_] <= '9') || s_[i_] == '.' || s_[i_] == 'e' || s_[i_] == 'E' ||
                                      s_[i_] == '+' || s_[i_] == '-')) ++i_;
            v.text = s_.substr(b, i_ - b);
        } else {
            fail("unexpected character");
        }
        return v;
    }

    const std::string& s_;
    std::size_t i_ = 0;
};

inline Value parse(const std::string& s) { return Parser(s).parse(); }

// A string as a JSON literal, quotes included.
inline std::string quote(const std::string& s) {
    std::string o = "\"";
    for (const char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof buf, "\\u%04x", static_cast<unsigned>(static_cast<unsigned char>(c)));
                    o += buf;
                } else {
                    o += c;
                }
        }
    }
    return o + "\"";
}

}  // namespace panvar::json
