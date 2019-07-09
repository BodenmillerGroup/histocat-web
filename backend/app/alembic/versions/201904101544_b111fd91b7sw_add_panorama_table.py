"""Add panorama table

Revision ID: b111fd91b7sw
Revises: a596fd31b7bd
Create Date: 2019-04-10 15:44:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.

revision = 'b111fd91b7sw'
down_revision = 'a596fd31b7bd'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'panorama',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('slide_id', sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True),
        sa.Column('metaname', sa.String(4096), index=True),
        sa.Column('original_id', sa.Integer(), index=True),
        sa.Column('description', sa.String()),
        sa.Column('slide_x4_pos_um', sa.Float()),
        sa.Column('slide_y4_pos_um', sa.Float()),
        sa.Column('slide_x3_pos_um', sa.Float()),
        sa.Column('slide_y3_pos_um', sa.Float()),
        sa.Column('slide_x2_pos_um', sa.Float()),
        sa.Column('slide_y2_pos_um', sa.Float()),
        sa.Column('slide_x1_pos_um', sa.Float()),
        sa.Column('slide_y1_pos_um', sa.Float()),
        sa.Column('image_end_offset', sa.BigInteger()),
        sa.Column('image_start_offset', sa.BigInteger()),
        sa.Column('pixel_width', sa.Integer()),
        sa.Column('pixel_height', sa.Integer()),
        sa.Column('image_format', sa.String(16)),
        sa.Column('pixel_scale_coef', sa.Float()),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('panorama')
