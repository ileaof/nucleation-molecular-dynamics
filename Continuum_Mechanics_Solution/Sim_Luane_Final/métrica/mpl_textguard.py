# -*- coding: utf-8 -*-
"""
Migration shim for matplotlib >= 3.10.

matplotlib 3.10 rewrote the ``ft2font`` C extension with pybind11, so
``FT2Font._set_transform(matrix, delta)`` now rejects a ``delta`` that does
not fit a C ``int32`` (it raises ``TypeError: _set_transform(): incompatible
function arguments``).  A ``Text`` artist placed far off the axes -- e.g. a
decorative alloy label drawn at a large negative data-y on a figure whose
y-axis is not scaled to match -- lands tens of millions of pixels outside the
figure.  freetype receives that pixel position as a 26.6 fixed-point delta
(``px * 64``), which overflows ``int32`` at about 3.36e7 px.  Older matplotlib
silently wrapped/ignored the overflow and simply did not draw the off-screen
text; the new binding makes it fatal, so the whole figure fails to render.

Importing this module wraps the Agg renderer's ``draw_text`` so that any text
whose display position would overflow ``int32`` is skipped -- exactly the old
behaviour.  The skipped text was always off-screen and invisible, so this
changes nothing the user sees; it only stops the crash.  No call site needs to
change.

    import mpl_textguard   # near the top, after importing matplotlib

The guard is idempotent and only affects text rendered > ~33 million px outside
the figure, so it is safe to import in every script (on-chart text is untouched).
"""

_INT32 = 2147483647
# freetype encodes the pixel translation as a 26.6 fixed-point delta (px * 64);
# that overflows int32 when |px| > INT32_MAX / 64 ~ 3.355e7.  Skip such text.
_PX_LIMIT = 3.3e7


def _install():
    try:
        from matplotlib.backends import backend_agg
    except Exception:
        return
    RendererAgg = backend_agg.RendererAgg
    if getattr(RendererAgg.draw_text, "_textguard_installed", False):
        return
    orig = RendererAgg.draw_text

    def safe_draw_text(self, gc, x, y, s, *args, **kwargs):
        try:
            if abs(float(x)) > _PX_LIMIT or abs(float(y)) > _PX_LIMIT:
                return  # off-screen text whose ft2font delta would overflow int32
        except Exception:
            pass
        try:
            return orig(self, gc, x, y, s, *args, **kwargs)
        except (TypeError, RuntimeError) as e:
            # mpl>=3.10 ft2font: int32 overflow (off-chart text) or freetype
            # raster overflow (a glyph bitmap too large to render). Both come
            # from text far outside the figure; old mpl skipped it silently.
            if any(k in str(e) for k in ("_set_transform", "FT_Render_Glyph", "raster overflow")):
                return
            raise

    safe_draw_text._textguard_installed = True
    RendererAgg.draw_text = safe_draw_text


_install()