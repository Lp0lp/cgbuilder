/* ===========================================================================
   DOM helper
   ===========================================================================
   Typed getElementById. The app's markup is fixed and known at author time,
   so the elements this looks up are treated as guaranteed-present and typed
   as their concrete element class (pass the element type as the generic, e.g.
   `byId<HTMLInputElement>('toggle-cg')`). Handlers that legitimately guard
   against a missing element (`if (el) ...`) still work at runtime — the guard
   simply reads as always-true to the type checker. */

/**
 * `document.getElementById(id)`, cast to the requested element type and
 * treated as non-null.
 * @param id - element id
 * @returns the element, typed as T
 */
export function byId<T extends HTMLElement = HTMLElement>(id: string): T {
  return document.getElementById(id) as T;
}
