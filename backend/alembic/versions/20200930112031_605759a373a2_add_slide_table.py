"""Add slide table

Revision ID: 605759a373a2
Revises: b33b3a0d8331
Create Date: 2020-09-30 11:20:31.970564

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = '605759a373a2'
down_revision = 'b33b3a0d8331'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'slide',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('project_id', sa.Integer(), sa.ForeignKey("project.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('origin_id', sa.Integer(), index=True),
        sa.Column('name', sa.String(4096), index=True),
        sa.Column('width_um', sa.Integer()),
        sa.Column('height_um', sa.Integer()),
        sa.Column('has_slide_image', sa.Boolean(), nullable=False, default=False),
        sa.Column('meta', JSONB()),
        sa.Column('session_meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )
    op.create_unique_constraint(
        'uq_slide_project_id_and_name',
        'slide',
        ['project_id', 'name']
    )


def downgrade():
    op.drop_table('slide')
    op.drop_unique_constraint('uq_slide_project_id_and_name')
