<template>
  <v-col>
    <v-toolbar dense class="toolbar">
      <v-toolbar-title>
        Group Members
      </v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn v-if="isGroupAdmin" text :to="`/main/groups/${activeGroupId}/members/create`" color="primary">
          Create Member
        </v-btn>
      </v-toolbar-items>
    </v-toolbar>
    <v-expansion-panels>
      <v-expansion-panel>
        <v-expansion-panel-header>Filter</v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-select
            v-model="roleFilter"
            :items="roles"
            item-text="text"
            item-value="value"
            chips
            clearable
            label="Role"
            multiple
            prepend-icon="mdi-filter-outline"
            solo
            dense
          >
            <template v-slot:selection="{ attrs, item, select, selected }">
              <v-chip
                v-bind="attrs"
                :input-value="selected"
                close
                @click="select"
                @click:close="removeRoleFilter(item)"
              >
                {{ item.text }}
              </v-chip>
            </template>
          </v-select>
          <v-switch label="Show inactive members" v-model="showInactiveMembers" />
        </v-expansion-panel-content>
      </v-expansion-panel>
    </v-expansion-panels>
    <v-card>
      <v-card-title>
        <v-spacer />
        <v-text-field v-model="search" append-icon="mdi-magnify" label="Search" single-line hide-details clearable />
      </v-card-title>
      <v-data-table
        :headers="headers"
        :items="items"
        :loading="!items"
        :search="search"
        :custom-filter="filter"
        :items-per-page="15"
        :footer-props="{
          itemsPerPageOptions: [10, 15, 20, -1],
          showFirstLastPage: true,
          showCurrentPage: true,
        }"
        multi-sort
      >
        <template v-slot:item.user="{ item }">
          <router-link
            v-if="item.user"
            class="link"
            :to="{
              name: 'main-admin-users-edit',
              params: {
                groupId: activeGroupId,
                id: item.user.id,
              },
            }"
          >
            {{ item.user.name }}
          </router-link>
        </template>
        <template v-slot:item.role="{ item }">
          {{ roleToString(item.role) }}
        </template>
        <template v-slot:item.isActive="{ item }">
          <v-icon v-if="item.isActive">mdi-check</v-icon>
        </template>
        <template v-slot:item.allPanels="{ item }">
          <v-icon v-if="item.allPanels">mdi-check</v-icon>
        </template>
        <template v-slot:item.action="{ item }">
          <v-menu bottom left>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on">
                <v-icon>mdi-dots-vertical</v-icon>
              </v-btn>
            </template>
            <v-list dense>
              <v-list-item
                :to="{
                  name: 'main-group-members-edit',
                  params: {
                    groupId: activeGroupId,
                    id: item.id,
                  },
                }"
              >
                <v-list-item-icon>
                  <v-icon color="grey">mdi-pencil-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Edit</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
              <v-list-item v-if="isGroupAdmin" @click="deleteMember(item.id)">
                <v-list-item-icon>
                  <v-icon color="red accent-1">mdi-delete-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Delete</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
            </v-list>
          </v-menu>
        </template>
      </v-data-table>
    </v-card>
  </v-col>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { memberModule } from "@/modules/member";
import { roleToString } from "@/utils/converters";
import { roleEnum } from "@/utils/enums";

@Component
export default class MembersListView extends Vue {
  readonly groupContext = groupModule.context(this.$store);
  readonly memberContext = memberModule.context(this.$store);

  readonly roleToString = roleToString;

  readonly roles = roleEnum;

  readonly headers = [
    // {
    //   text: "Id",
    //   sortable: true,
    //   value: "id",
    //   align: "right",
    //   filterable: false,
    //   width: "80",
    // },
    {
      text: "Name",
      value: "user",
      sort: (a, b) => {
        if (a === null) {
          return 1;
        }
        if (b === null) {
          return -1;
        }
        return a.name.localeCompare(b.name);
      },
    },
    {
      text: "Role",
      value: "role",
    },
    {
      text: "Active",
      sortable: true,
      value: "isActive",
      align: "left",
      filterable: false,
      width: "100",
    },
    {
      text: "All panels",
      sortable: true,
      value: "allPanels",
      align: "left",
      filterable: false,
      width: "130",
    },
    {
      text: "Actions",
      value: "action",
      sortable: false,
      filterable: false,
      width: "70",
    },
  ];

  search = "";

  roleFilter: number[] = [];
  showInactiveMembers = false;

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  get isGroupAdmin() {
    return this.groupContext.getters.isGroupAdmin;
  }

  get items() {
    let items = this.memberContext.getters.members;
    if (!this.showInactiveMembers) {
      items = items.filter((item) => item.is_active);
    }
    if (this.roleFilter.length > 0) {
      items = items.filter((item) => this.roleFilter.includes(item.role));
    }
    return items;
  }

  filter(value, search, item) {
    if (!search) {
      return true;
    }
    const normalizedSearchTerm = search.toLowerCase().trim();
    return item.user ? item.user.name.toLowerCase().indexOf(normalizedSearchTerm) !== -1 : false;
  }

  async mounted() {
    await this.memberContext.actions.getGroupMembers(+this.$router.currentRoute.params.groupId);
  }

  async deleteMember(id: number) {
    if (self.confirm("Are you sure you want to delete the group member?")) {
      await this.memberContext.actions.deleteMember(id);
    }
  }

  removeRoleFilter(item) {
    this.roleFilter.splice(this.roleFilter.indexOf(item), 1);
    this.roleFilter = [...this.roleFilter];
  }
}
</script>

<style scoped>
.toolbar {
  margin-bottom: 10px;
}
</style>
