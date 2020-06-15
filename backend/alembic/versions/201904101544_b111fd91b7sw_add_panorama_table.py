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
        sa.Column('origin_id', sa.Integer(), index=True),
        sa.Column('image_type', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('x1', sa.Float()),
        sa.Column('y1', sa.Float()),
        sa.Column('x2', sa.Float()),
        sa.Column('y2', sa.Float()),
        sa.Column('x3', sa.Float()),
        sa.Column('y3', sa.Float()),
        sa.Column('x4', sa.Float()),
        sa.Column('y4', sa.Float()),
        sa.Column('rotation_angle', sa.Float()),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
    )


def downgrade():
    op.drop_table('panorama')
