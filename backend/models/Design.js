const mongoose = require('mongoose');

const DesignSchema = new mongoose.Schema({
  name:        { type: String, default: 'Untitled Wing' },
  mode:        { type: String, enum: ['parametric','freeform'], default: 'parametric' },

  /* parametric wing fields */
  span:        { type: Number },
  root_chord:  { type: Number },
  tip_chord:   { type: Number },
  sweep_deg:   { type: Number },
  dihedral_deg:{ type: Number },

  /* freeform polygon fields */
  points:      [{ x: Number, y: Number }],

  /* aerodynamic results at save time */
  CL:          { type: Number },
  CD:          { type: Number },
  LD_ratio:    { type: Number },

  /* meta */
  author:      { type: String, default: 'Anonymous' },
  createdAt:   { type: Date,   default: Date.now }
});

module.exports = mongoose.model('Design', DesignSchema);