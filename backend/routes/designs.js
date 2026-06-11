const express  = require('express');
const router   = express.Router();
const Design   = require('../models/Design');

/* ── Save a design ── POST /api/designs */
router.post('/', async (req, res) => {
  try {
    const design = new Design(req.body);
    await design.save();
    res.json({ success: true, id: design._id, design });
  } catch (err) {
    res.status(500).json({ success: false, error: err.message });
  }
});

/* ── Load a design ── GET /api/designs/:id */
router.get('/:id', async (req, res) => {
  try {
    const design = await Design.findById(req.params.id);
    if (!design) return res.status(404).json({ success: false, error: 'Not found' });
    res.json({ success: true, design });
  } catch (err) {
    res.status(500).json({ success: false, error: err.message });
  }
});

/* ── Get all designs ── GET /api/designs */
router.get('/', async (req, res) => {
  try {
    const designs = await Design.find()
      .sort({ createdAt: -1 })
      .limit(20)
      .select('name mode CL LD_ratio author createdAt');
    res.json({ success: true, designs });
  } catch (err) {
    res.status(500).json({ success: false, error: err.message });
  }
});

/* ── Delete a design ── DELETE /api/designs/:id */
router.delete('/:id', async (req, res) => {
  try {
    await Design.findByIdAndDelete(req.params.id);
    res.json({ success: true });
  } catch (err) {
    res.status(500).json({ success: false, error: err.message });
  }
});

module.exports = router;