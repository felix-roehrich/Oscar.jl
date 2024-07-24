struct SparseWeylGroup
  root_system::RootSystem
  cartan_matrix::SMat
  cartan_matrix_t::SMat
end

struct SparseWeylGroupElem
  # parent group
  parent::WeylGroup

  # short revlex normal form of the word
  word::Vector{UInt8}
end

function Base.:(*)(x::SparseWeylGroupElem, y::SparseWeylGroupElem)
  @req x.parent === y.parent "$x, $y must belong to the same Weyl group"

  word = copy(y.word)
  for s in Iterators.reverse(x.word)
    lmul!(x.parent.refl, word, s)
  end

  return WeylGroupElem(x.parent, word)
end

function sparse_lmul!(R::AbstractRootSystem, cartan_matrix::SMat, word::Vector{Int}, s::Int)
  insert_index = 1
  insert_letter = s

  root = r.root_system.positive_roots[s]
  for i in 1:length(word)
    if word[i] == root
      deleteat!(word, i)
      return word
    end

    word[i], root
    addmul!(root, r.root_system.positive_roots[i], dot(cartan_matrix[s], r.vec))
    root = refl[Int(word[i]), Int(root)]
    if root == 0
      # r is no longer a minimal root, meaning we found the best insertion point
      break
    end

    # check if we have a better insertion point now. Since word[i] is a simple
    # root, if root < word[i] it must be simple.
    if root < word[i]
      insert_index = i + 1
      insert_letter = T(root)
    end
  end

  insert!(word, insert_index, insert_letter)
  return
end
