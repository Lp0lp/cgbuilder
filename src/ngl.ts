/* ===========================================================================
   NGL adapter
   ===========================================================================
   Re-exposes the pieces of the `ngl` npm package the app uses under the
   minimal `NGLStatic` surface it already relies on, replacing the former
   `NGL` browser global (which was loaded via a <script> tag). NGL re-exports
   three's `Vector3`, so all three constructors come from the one package.

   All the type impedance between our deliberately-partial interfaces (see
   types.ts) and NGL's full, much larger types is isolated to the single cast
   here — consumers keep using `NGL.Stage` / `NGL.Shape` / `NGL.Vector3`
   exactly as before. */
import { Stage, Shape, Vector3 } from 'ngl';
import type { NGLStatic } from './types.js';

export const NGL = { Stage, Shape, Vector3 } as unknown as NGLStatic;
