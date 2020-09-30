"""Add panorama table

Revision ID: 69e599da5845
Revises: 605759a373a2
Create Date: 2020-09-30 11:21:17.093851

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = '69e599da5845'
down_revision = '605759a373a2'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'panorama',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('slide_id', sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True, nullable=False),
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
