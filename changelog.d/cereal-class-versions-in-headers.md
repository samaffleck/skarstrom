### Fixed — class versions visible wherever a deck is saved or loaded
- Every serialised engine class now declares its `CEREAL_CLASS_VERSION` in the
  header that holds its inline `serialize` template, not in its `.cpp`. A
  translation unit that saved or loaded a deck without seeing the declaration
  used version 0, so which version a save wrote (and which branch a load took)
  depended on which instantiation the linker kept. Saved decks do not change.
  A new unit test (`CerealClassVersions.*`) pins the version each class writes.
