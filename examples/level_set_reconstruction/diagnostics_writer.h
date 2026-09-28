#ifndef DIAGNOSTICS_WRITER_H
#define DIAGNOSTICS_WRITER_H

#include <fstream>
#include <functional>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

class DiagnosticsWriter {
 public:
  explicit DiagnosticsWriter(const std::string& filename) : out(filename) {
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out << std::setprecision(17);
  }

  template <class T>
  void add(const std::string& name, const T& value) {
    columns.push_back({name, [&value](std::ostream& out) { out << value; }});
  }

  void writeHeader() {
    for (std::size_t i = 0; i < columns.size(); ++i) {
      if (i > 0) out << ' ';
      out << columns[i].name;
    }
    out << '\n';
  }

  void write() {
    for (std::size_t i = 0; i < columns.size(); ++i) {
      if (i > 0) out << ' ';
      columns[i].write(out);
    }
    out << '\n';
  }

 private:
  struct Column {
    std::string name;
    std::function<void(std::ostream&)> write;
  };

  std::vector<Column> columns;
  std::ofstream out;
};

#endif