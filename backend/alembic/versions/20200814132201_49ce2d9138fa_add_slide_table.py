"""Add slide table

Revision ID: 49ce2d9138fa
Revises: af7cd5a143a9
Create Date: 2020-08-14 13:22:01.723775

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = '49ce2d9138fa'
down_revision = 'af7cd5a143a9'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'slide',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('project_id', sa.Integer(), sa.ForeignKey("project.id", ondelete="CASCADE"), index=True),
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

