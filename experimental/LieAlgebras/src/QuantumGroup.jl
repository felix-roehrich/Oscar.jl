struct QuantumGroup
  root_system::RootSystem
end

function quantum_group(R::RootSystem)
  return QuantumGroup(R)
end
