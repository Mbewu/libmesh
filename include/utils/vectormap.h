// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef LIBMESH_VECTORMAP_H
#define LIBMESH_VECTORMAP_H

// C++ Includes   -----------------------------------
#include <vector>



namespace libMesh
{

/**
 * This \p vectormap templated class is intended to provide the
 * performance characteristics of a sorted std::vector with an
 * interface more closely resembling that of a std::map, for use in
 * particular when memory is tight.
 *
 * \author  Benjamin S. Kirk
 */



  template <typename Key, typename Tp,
	    typename Compare = std::less<Key>,
	    typename Alloc = std::allocator<std::pair<const Key, Tp> > >
  class vectormap : public std::vector<std::pair<Key, Tp>, Alloc >
  {
  public:

    typedef Key                      key_type;
    typedef Tp                       mapped_type;
    typedef std::pair<const Key, Tp> value_type;
    typedef Compare                  key_compare;
    typedef Alloc                    allocator_type;

    /**
     * Default constructor.  Initializes sorted member to false.
     */
    vectormap () :
      _is_sorted(false)
    {}

    /**
     * Inserts \p x into the vectormap.
     */
    void insert (const value_type &x)
    {
      _is_sorted = false;
      this->push_back(x);
    }

  private:

    bool _is_sorted;
  };


// template <typename Val, typename index_t=unsigned int>
// class mapvector : public std::map<index_t, Val>
// {
// public:
//   typedef std::map<index_t, Val> maptype;

//   Val& operator[] (const index_t &k)
//   {
//     return maptype::operator[](k);
//   }
//   Val operator[] (const index_t &k) const
//   {
//     typename maptype::const_iterator it = this->find(k);
//       return it == this->end().it? Val() : it->second;
//   }

//   class veclike_iterator
//   {
//   public:
//     veclike_iterator(const typename maptype::iterator &i)
//       : it(i) {}

//     veclike_iterator(const veclike_iterator &i)
//       : it(i.it) {}

//     Val& operator*() const { return it->second; }

//     veclike_iterator& operator++() { ++it; return *this; }

//     veclike_iterator operator++(int) {
//       veclike_iterator i = *this;
//       ++(*this);
//       return i;
//     }

//     bool operator==(const veclike_iterator &other) const {
//       return it == other.it;
//     }

//     bool operator!=(const veclike_iterator &other) const {
//       return it != other.it;
//     }

//     typename maptype::iterator it;
//   };

//   class const_veclike_iterator
//   {
//   public:
//     const_veclike_iterator(const typename maptype::const_iterator &i)
//       : it(i) {}

//     const_veclike_iterator(const const_veclike_iterator &i)
//       : it(i.it) {}

//     const_veclike_iterator(const veclike_iterator &i)
//       : it(i.it) {}

//     const Val& operator*() const { return it->second; }

//     const_veclike_iterator& operator++() { ++it; return *this; }

//     const_veclike_iterator operator++(int) {
//       veclike_iterator i = *this;
//       ++(*this);
//       return i;
//     }

//     bool operator==(const const_veclike_iterator &other) const {
//       return it == other.it;
//     }

//     bool operator!=(const const_veclike_iterator &other) const {
//       return it != other.it;
//     }

//     typename maptype::const_iterator it;
//   };

//   void erase(index_t i) {
//       maptype::erase(i);
//   }

//   void erase(const veclike_iterator &pos) {
//       maptype::erase(pos.it);
//   }

//   veclike_iterator begin() {
//     return veclike_iterator(maptype::begin());
//   }

//   const_veclike_iterator begin() const {
//     return const_veclike_iterator(maptype::begin());
//   }

//   veclike_iterator end() {
//     return veclike_iterator(maptype::end());
//   }

//   const_veclike_iterator end() const {
//     return const_veclike_iterator(maptype::end());
//   }
// };

} // namespace libMesh

#endif // LIBMESH_VECTORMAP_H
